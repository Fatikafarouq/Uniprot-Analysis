from __future__ import annotations

import re
from typing import Any

from .evidence import evidence_label


def _clean(value: Any) -> str:
    return " ".join(str(value or "").split()).strip()


def _ensure_period(text: str) -> str:
    text = _clean(text)
    if not text:
        return text
    return text if text[-1] in ".!?" else text + "."


def _split_sentences(text: str) -> list[str]:
    text = _clean(text)
    if not text:
        return []
    parts = re.split(r"(?<=[.!?])\s+(?=[A-Z0-9])", text)
    return [_ensure_period(part) for part in parts if _clean(part)]


def _source(field: str, text: str, evidences: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    return {
        "field": field,
        "text": _clean(text),
        "evidence": [evidence_label(item) for item in (evidences or []) if isinstance(item, dict)],
    }


def _comment_texts(record: dict[str, Any], comment_type: str) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for comment in record.get("comments", []) or []:
        if comment.get("commentType") != comment_type:
            continue
        for item in comment.get("texts", []) or []:
            if isinstance(item, dict) and item.get("value"):
                output.append(_source(comment_type, str(item["value"]), item.get("evidences", []) or []))
    return output


def _subcellular_locations(record: dict[str, Any]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    seen: set[str] = set()
    for comment in record.get("comments", []) or []:
        if comment.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        for block in comment.get("subcellularLocations", []) or []:
            if not isinstance(block, dict):
                continue
            location = block.get("location") or {}
            topology = block.get("topology") or {}
            orientation = block.get("orientation") or {}
            values = []
            evidences: list[dict[str, Any]] = []
            for node in (location, topology, orientation):
                if isinstance(node, dict) and node.get("value"):
                    values.append(str(node["value"]))
                    evidences.extend(node.get("evidences", []) or [])
            text = "; ".join(values)
            key = text.lower()
            if text and key not in seen:
                seen.add(key)
                output.append(_source("SUBCELLULAR LOCATION", text, evidences))
    return output


def _catalytic_activities(record: dict[str, Any]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for comment in record.get("comments", []) or []:
        if comment.get("commentType") != "CATALYTIC ACTIVITY":
            continue
        reaction = comment.get("reaction") or {}
        if not isinstance(reaction, dict):
            continue
        name = reaction.get("name")
        if name:
            output.append(_source("CATALYTIC ACTIVITY", str(name), reaction.get("evidences", []) or []))
    return output


def _cofactors(record: dict[str, Any]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for comment in record.get("comments", []) or []:
        if comment.get("commentType") != "COFACTOR":
            continue
        for cofactor in comment.get("cofactors", []) or []:
            if not isinstance(cofactor, dict):
                continue
            name = cofactor.get("name")
            description = cofactor.get("description")
            if name:
                text = str(name)
                if description:
                    text += f": {description}"
                output.append(_source("COFACTOR", text, cofactor.get("evidences", []) or []))
    return output


def _binding_features(record: dict[str, Any]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for feature in record.get("features", []) or []:
        feature_type = str(feature.get("type") or "")
        if feature_type not in {"Binding site", "Active site", "Site", "Metal binding"}:
            continue
        description = _clean(feature.get("description"))
        if not description:
            continue
        output.append(_source(f"{feature_type.upper()} FEATURE", description, feature.get("evidences", []) or []))
    return output


def _friendly_location(text: str) -> str:
    lower = text.lower()
    if "cell membrane" in lower or "plasma membrane" in lower:
        return "the cell membrane"
    if "nucleus" in lower:
        return "the nucleus"
    if "cytoplasm" in lower or "cytosol" in lower:
        return "the cytoplasm"
    if "mitochond" in lower:
        return "the mitochondria"
    if "endoplasmic reticulum" in lower:
        return "the endoplasmic reticulum"
    if "golgi" in lower:
        return "the Golgi apparatus"
    if "secreted" in lower or "extracellular" in lower:
        return "outside the cell"
    return text


def _friendly_function_sentence(text: str) -> str:
    sentence = _ensure_period(text)
    bare = sentence[:-1] if sentence.endswith(".") else sentence

    # Common UniProt receptor wording. This is a wording transformation only:
    # every biological relationship remains present in the source sentence.
    match = re.match(
        r"^Receptor tyrosine kinase which mediates the (?:pleiotropic )?actions of (.+)$",
        bare,
        flags=re.I,
    )
    if match:
        ligand = _clean(match.group(1))
        return f"This protein is a receptor tyrosine kinase that helps cells respond to {ligand}."

    match = re.match(r"^Receptor for (.+)$", bare, flags=re.I)
    if match:
        return f"This protein acts as a receptor for {_clean(match.group(1))}."

    match = re.match(r"^Binding of (.+?) leads to (.+)$", bare, flags=re.I)
    if match:
        ligand = _clean(match.group(1))
        consequence = _clean(match.group(2))
        return f"When {ligand} binds, {consequence[0].lower() + consequence[1:] if consequence else consequence}."

    match = re.match(r"^Upon binding of (.+?),? (.+)$", bare, flags=re.I)
    if match:
        ligand = _clean(match.group(1))
        consequence = _clean(match.group(2))
        return f"When {ligand} binds, {consequence[0].lower() + consequence[1:] if consequence else consequence}."

    replacements = (
        (r"^Catalyzes (.+)$", "This protein is an enzyme that catalyzes {x}."),
        (r"^Transports (.+)$", "This protein transports {x}."),
        (r"^Binds (.+)$", "This protein binds {x}."),
        (r"^Required for (.+)$", "This protein is needed for {x}."),
        (r"^Involved in (.+)$", "This protein helps with {x}."),
        (r"^Acts as (.+)$", "This protein acts as {x}."),
        (r"^Functions as (.+)$", "This protein functions as {x}."),
        (r"^Mediates (.+)$", "This protein helps mediate {x}."),
    )
    for pattern, template in replacements:
        match = re.match(pattern, bare, flags=re.I)
        if match:
            return template.format(x=_clean(match.group(1)))

    # If the source already gives a complete mechanistic sentence, preserve it
    # rather than inventing a simpler relationship that UniProt did not state.
    return sentence


def _friendly_regulation(text: str) -> str:
    sentence = _ensure_period(text)
    bare = sentence[:-1] if sentence.endswith(".") else sentence

    match = re.match(r"^Activated by (.+)$", bare, flags=re.I)
    if match:
        return f"It is activated by {_clean(match.group(1))}."
    match = re.match(r"^Inhibited by (.+)$", bare, flags=re.I)
    if match:
        return f"It is inhibited by {_clean(match.group(1))}."
    match = re.match(r"^Binding of (.+?) leads to (.+)$", bare, flags=re.I)
    if match:
        consequence = _clean(match.group(2))
        return f"When {_clean(match.group(1))} binds, {consequence[0].lower() + consequence[1:] if consequence else consequence}."
    return sentence


def _combine_location_with_role(location: str, role: str) -> str | None:
    if location != "the cell membrane":
        return None
    if role.startswith("This protein is a receptor"):
        return role.replace("This protein is", "This protein sits in the cell membrane and is", 1)
    if role.startswith("This protein acts as a receptor"):
        return role.replace("This protein acts", "This protein sits in the cell membrane and acts", 1)
    return None


def mechanism_explanation(record: dict[str, Any]) -> dict[str, Any]:
    """Build a beginner-facing mechanism explanation using only UniProt fields.

    This function is intentionally conservative. It can rephrase or combine
    relationships that UniProt states, but it never supplies a missing causal
    step from general biological knowledge.
    """
    function_items = _comment_texts(record, "FUNCTION")
    regulation_items = _comment_texts(record, "ACTIVITY REGULATION")
    subunit_items = _comment_texts(record, "SUBUNIT")
    pathway_items = _comment_texts(record, "PATHWAY")
    locations = _subcellular_locations(record)
    catalytic = _catalytic_activities(record)
    cofactors = _cofactors(record)
    binding = _binding_features(record)

    summary: list[str] = []
    used_sources: list[dict[str, Any]] = []

    function_sentences: list[tuple[str, dict[str, Any]]] = []
    for item in function_items:
        for sentence in _split_sentences(item["text"]):
            function_sentences.append((sentence, item))

    role = None
    role_source = None
    if function_sentences:
        role = _friendly_function_sentence(function_sentences[0][0])
        role_source = function_sentences[0][1]

    friendly_location = _friendly_location(locations[0]["text"]) if locations else None
    combined = _combine_location_with_role(friendly_location, role) if (friendly_location and role) else None

    if combined:
        summary.append(combined)
        used_sources.extend([locations[0], role_source])
    elif role:
        summary.append(role)
        used_sources.append(role_source)
        if friendly_location:
            if friendly_location == "outside the cell":
                summary.append("It is secreted outside the cell.")
            elif friendly_location.startswith("the "):
                summary.append(f"It is found in {friendly_location}.")
            else:
                summary.append(f"UniProt places it at {friendly_location}.")
            used_sources.append(locations[0])
    elif friendly_location:
        if friendly_location == "outside the cell":
            summary.append("This protein is secreted outside the cell.")
        elif friendly_location.startswith("the "):
            summary.append(f"This protein is found in {friendly_location}.")
        else:
            summary.append(f"UniProt places this protein at {friendly_location}.")
        used_sources.append(locations[0])

    # A second FUNCTION sentence often contains the actual trigger/consequence
    # (for example ligand binding -> phosphorylation -> signalling).
    for sentence, source in function_sentences[1:3]:
        friendly = _friendly_function_sentence(sentence)
        if friendly and friendly not in summary:
            summary.append(friendly)
            used_sources.append(source)
        if len(summary) >= 3:
            break

    if len(summary) < 3 and regulation_items:
        friendly = _friendly_regulation(regulation_items[0]["text"])
        if friendly and friendly not in summary:
            summary.append(friendly)
            used_sources.append(regulation_items[0])

    if len(summary) < 3 and catalytic:
        summary.append(f"It catalyzes this reaction: {catalytic[0]['text']}.")
        used_sources.append(catalytic[0])

    if len(summary) < 4 and pathway_items:
        summary.append(f"It participates in this pathway: {_ensure_period(pathway_items[0]['text'])}")
        used_sources.append(pathway_items[0])

    # Add an interaction only when it adds a mechanistic relation that has not
    # already been captured by the FUNCTION/REGULATION wording.
    if len(summary) < 4 and subunit_items:
        interaction = _ensure_period(subunit_items[0]["text"])
        lower = interaction.lower()
        if any(word in lower for word in ("bind", "interact", "complex", "dimer", "subunit")):
            if not any(_clean(interaction).lower() in _clean(item).lower() for item in summary):
                summary.append(interaction)
                used_sources.append(subunit_items[0])

    if len(summary) < 4 and cofactors:
        name = cofactors[0]["text"]
        summary.append(f"UniProt reports this cofactor requirement: {name}.")
        used_sources.append(cofactors[0])

    if len(summary) < 4 and binding:
        summary.append(f"UniProt annotates a functional site as: {binding[0]['text']}.")
        used_sources.append(binding[0])

    # Deduplicate sources by field + exact source wording.
    deduped_sources: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for item in used_sources:
        if not item:
            continue
        key = (str(item.get("field")), str(item.get("text")))
        if key in seen:
            continue
        seen.add(key)
        deduped_sources.append(item)

    evidence_categories = {
        evidence.get("category")
        for item in deduped_sources
        for evidence in item.get("evidence", []) or []
        if isinstance(evidence, dict)
    }
    if "experimental" in evidence_categories:
        evidence_note = "At least one UniProt annotation used here includes experimental evidence."
    elif "similarity" in evidence_categories:
        evidence_note = "The UniProt annotations used here include similarity-based evidence."
    elif "computational" in evidence_categories:
        evidence_note = "The UniProt annotations used here include computational or sequence-model evidence."
    elif deduped_sources:
        evidence_note = "This explanation is a plain-language synthesis of the UniProt annotations listed below."
    else:
        evidence_note = "UniProt currently provides limited mechanistic detail for this entry."

    if not summary:
        summary = ["UniProt currently provides limited mechanistic detail for this entry."]

    return {
        "summary": summary[:4],
        "sources": deduped_sources,
        "evidence_note": evidence_note,
        "limited": not bool(deduped_sources),
    }
