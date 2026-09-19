from __future__ import annotations

import re
from typing import Any

from .records import get_all_gene_labels, get_all_protein_search_names, get_alternative_names, meaningful_query_tokens
from .text import normalize_text

SEMANTIC_XREF_KEY_WORDS = {"term", "pathway", "disease", "phenotype", "process", "function"}


def literature_title_identifies_record(record: dict[str, Any], title: str) -> bool:
    """A literature title counts only if it identifies this record's protein/gene."""
    title_norm = normalize_text(title)
    if not title_norm:
        return False

    labels = [*get_all_gene_labels(record), *get_all_protein_search_names(record), *get_alternative_names(record)]
    seen: set[str] = set()
    for label in labels:
        label_norm = normalize_text(label)
        if not label_norm or label_norm in seen or len(label_norm) < 3:
            continue
        seen.add(label_norm)
        pattern = r"(?<![a-z0-9])" + re.escape(label_norm) + r"(?![a-z0-9])"
        if re.search(pattern, title_norm):
            return True
    return False


def recursive_strings(value: Any):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for child in value.values():
            yield from recursive_strings(child)
    elif isinstance(value, list):
        for child in value:
            yield from recursive_strings(child)


def iter_discovery_strings(record: dict[str, Any]):
    """Port of the Colab v14 discovery evidence surface."""
    for name in get_all_protein_search_names(record):
        if name:
            yield "Protein name", str(name)

    for label in get_all_gene_labels(record):
        if label:
            yield "Gene name or synonym", str(label)

    for comment in record.get("comments", []) or []:
        comment_type = str(comment.get("commentType", "ANNOTATION")).replace("_", " ")
        for value in recursive_strings(comment):
            if value:
                yield comment_type, str(value)

    for keyword in record.get("keywords", []) or []:
        value = keyword.get("name") if isinstance(keyword, dict) else keyword
        if value:
            yield "UniProt keyword", str(value)

    for feature in record.get("features", []) or []:
        feature_type = str(feature.get("type", "Feature"))
        description = feature.get("description")
        if description:
            yield f"{feature_type} annotation", str(description)

    # Cross-reference values are eligible only when the property name itself
    # describes semantic biological evidence. Generic identifiers are ignored.
    for xref in record.get("uniProtKBCrossReferences", []) or []:
        database = str(xref.get("database", "database"))
        for prop in xref.get("properties", []) or []:
            if not isinstance(prop, dict):
                continue
            value = prop.get("value")
            key = prop.get("key")
            if not value or not key:
                continue
            normalized_key = normalize_text(key).replace(" ", "")
            if not any(word in normalized_key for word in SEMANTIC_XREF_KEY_WORDS):
                continue
            yield f"{database} cross-reference ({key})", str(value)

    for reference in record.get("references", []) or []:
        citation = reference.get("citation", {}) or {}
        title = citation.get("title")
        if title and literature_title_identifies_record(record, str(title)):
            yield "Reference title linked by UniProt", str(title)


def text_matches_query(text: str, query: str) -> bool:
    candidate = normalize_text(text)
    phrase = normalize_text(query)
    if not candidate or not phrase:
        return False
    if phrase in candidate:
        return True
    tokens = meaningful_query_tokens(phrase) or phrase.split()
    candidate_tokens = set(candidate.split())
    return bool(tokens) and all(token in candidate_tokens for token in tokens)


def shorten_match_text(text: str, phrase: str, limit: int = 260) -> str:
    text = " ".join(str(text).split())
    if len(text) <= limit:
        return text

    normalized_text = normalize_text(text)
    query = normalize_text(phrase)
    position = normalized_text.find(query)
    if position < 0:
        position = -1
        for token in meaningful_query_tokens(query):
            position = normalized_text.find(token)
            if position >= 0:
                break
    if position < 0:
        return text[: limit - 1] + "…"

    start = max(0, position - limit // 3)
    end = min(len(text), start + limit)
    snippet = text[start:end]
    if start > 0:
        snippet = "…" + snippet
    if end < len(text):
        snippet += "…"
    return snippet


def annotation_match_contexts(record: dict[str, Any], query: str, max_items: int = 3) -> list[dict[str, str]]:
    """One source field must support the query; evidence is never stitched across fields."""
    contexts: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for source, value in iter_discovery_strings(record):
        clean = " ".join(str(value).split())
        signature = (source, clean)
        if not clean or signature in seen:
            continue
        seen.add(signature)
        if not text_matches_query(clean, query):
            continue
        contexts.append({"source": source, "text": shorten_match_text(clean, query)})
        if len(contexts) >= max_items:
            break
    return contexts


def annotation_match_contexts_for_phrases(
    record: dict[str, Any], phrases: list[str] | tuple[str, ...] | str, max_items: int = 3
) -> tuple[list[dict[str, str]], str | None]:
    if isinstance(phrases, str):
        phrases = [phrases]
    for phrase in phrases:
        contexts = annotation_match_contexts(record, phrase, max_items=max_items)
        if contexts:
            return contexts, phrase
    return [], None


DISCOVERY_EVIDENCE_TYPES = {
    "name_field": {"label": "Direct protein / gene name match"},
    "source_organism": {"label": "Source organism match"},
    "disease": {"label": "UniProt disease annotation"},
    "function": {"label": "UniProt function annotation"},
    "allergen": {"label": "UniProt allergen annotation"},
    "pathway": {"label": "UniProt pathway annotation"},
    "literature": {"label": "Literature linked by UniProt"},
    "linked_database": {"label": "Linked database evidence"},
    "name_field": {"label": "UniProt protein / gene name evidence"},
    "other_uniprot": {"label": "Other UniProt annotation"},
}


def provenance_type_for_source(source: str | None) -> str:
    source = str(source or "").strip()
    upper = source.upper()
    if upper == "DISEASE":
        return "disease"
    if upper == "FUNCTION":
        return "function"
    if upper == "ALLERGEN":
        return "allergen"
    if upper == "PATHWAY":
        return "pathway"
    if source in {"Source organism scientific name", "Source organism common name"}:
        return "source_organism"
    if source == "Reference title linked by UniProt":
        return "literature"
    if " cross-reference" in source:
        return "linked_database"
    if source in {"Protein name", "Gene name or synonym"}:
        return "name_field"
    return "other_uniprot"


def classify_discovery_evidence(contexts: list[dict[str, str]]) -> str:
    if not contexts:
        return "other_uniprot"
    return provenance_type_for_source(contexts[0].get("source"))



def beginner_connection_explanation(evidence_type: str, phrase: str | None) -> str:
    phrase = str(phrase or "the search concept")
    explanations = {
        "name_field": (
            f'The protein or gene name itself contains "{phrase}". This is a name-level match.'
        ),
        "source_organism": (
            f'The UniProt source-organism name contains "{phrase}". This tells you where the protein comes from; '
            f'it does not by itself mean UniProt says the protein function is "{phrase}".'
        ),
        "disease": (
            f'A UniProt DISEASE annotation for this protein mentions "{phrase}". The connection comes from disease annotation evidence.'
        ),
        "function": (
            f"A UniProt FUNCTION annotation for this protein mentions \"{phrase}\". The connection comes from UniProt's description of the protein function."
        ),
        "allergen": (
            f'A UniProt ALLERGEN annotation for this protein mentions "{phrase}".'
        ),
        "pathway": (
            f'A UniProt PATHWAY annotation for this protein mentions "{phrase}".'
        ),
        "literature": (
            f'A publication title linked by UniProt mentions "{phrase}" and also identifies this protein or gene.'
        ),
        "linked_database": (
            f'A meaningful field in a database cross-reference carried by the UniProt record mentions "{phrase}".'
        ),
        "other_uniprot": (
            f'Another traceable UniProt annotation contains "{phrase}". The exact annotation type is shown below.'
        ),
    }
    return explanations.get(evidence_type, explanations["other_uniprot"])


def evidence_label(evidence: dict[str, Any]) -> dict[str, Any]:
    code = str(evidence.get("evidenceCode") or "")
    source = evidence.get("source")
    source_id = evidence.get("id")

    if code == "ECO:0000269":
        category = "experimental"
        label = "Experimental evidence"
    elif code == "ECO:0000250":
        category = "similarity"
        label = "Inferred from sequence similarity"
    elif code in {"ECO:0000255", "ECO:0000256"}:
        category = "computational"
        label = "Sequence-model / computational evidence"
    elif code in {"ECO:0000312", "ECO:0000313"}:
        category = "imported"
        label = "Imported database evidence"
    else:
        category = "other"
        label = f"UniProt evidence code {code}" if code else "Evidence code not exposed"

    return {
        "code": code or None,
        "category": category,
        "label": label,
        "source": source,
        "source_id": source_id,
    }


def function_evidence(record: dict[str, Any]) -> dict[str, Any]:
    statements: list[dict[str, Any]] = []
    category_counts: dict[str, int] = {}

    for comment in record.get("comments", []) or []:
        if comment.get("commentType") != "FUNCTION":
            continue
        for text in comment.get("texts", []) or []:
            if not isinstance(text, dict) or not text.get("value"):
                continue
            evidences = [evidence_label(item) for item in text.get("evidences", []) or []]
            for item in evidences:
                category_counts[item["category"]] = category_counts.get(item["category"], 0) + 1
            statements.append({"text": text["value"], "evidence": evidences})

    if category_counts.get("experimental"):
        summary = "At least one UniProt function annotation has experimental evidence."
    elif category_counts.get("similarity"):
        summary = "The exposed function annotation is supported by similarity-based evidence rather than direct experimental evidence."
    elif category_counts.get("computational"):
        summary = "The exposed function annotation is supported by computational or sequence-model evidence."
    elif statements:
        summary = "UniProt exposes function annotation, but this prototype does not map its evidence codes to a stronger plain-language category."
    else:
        summary = "No UniProt FUNCTION comment was exposed for this record."

    return {"summary": summary, "category_counts": category_counts, "statements": statements}
