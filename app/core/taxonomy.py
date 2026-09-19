from __future__ import annotations

import re
from typing import Any

from .http import DataClient, SourceUnavailable
from .records import direct_match_details
from .text import normalize_text, singularize_phrase

SPECIFIC_RANKS = {
    "species",
    "subspecies",
    "strain",
    "isolate",
    "serotype",
    "forma",
    "varietas",
}

PLACEHOLDER_WORDS = {
    "uncultured",
    "unidentified",
    "unclassified",
    "environmental sample",
    "metagenome",
    "metagenomic",
}


def taxonomy_names(record: dict[str, Any]) -> list[tuple[str, str]]:
    output: list[tuple[str, str]] = []
    for label, key in (("scientific name", "scientificName"), ("common name", "commonName")):
        value = record.get(key)
        if value:
            output.append((label, str(value)))

    for key in ("synonyms", "otherNames"):
        for item in record.get(key, []) or []:
            if isinstance(item, str):
                output.append(("synonym", item))
            elif isinstance(item, dict):
                value = item.get("value") or item.get("name")
                if value:
                    output.append(("synonym", str(value)))
    return output


def taxonomy_display_name(record: dict[str, Any]) -> str:
    scientific = record.get("scientificName")
    common = record.get("commonName")
    if common and scientific and normalize_text(common) != normalize_text(scientific):
        return f"{common} ({scientific})"
    return str(common or scientific or record.get("taxonId") or "Unknown organism")


def taxonomy_record_is_placeholder(record: dict[str, Any]) -> bool:
    haystack = " ".join(value for _, value in taxonomy_names(record)).lower()
    return any(word in haystack for word in PLACEHOLDER_WORDS)


def _lookup_variants(phrase: str) -> list[str]:
    output: list[str] = []
    for value in (normalize_text(phrase), singularize_phrase(phrase)):
        value = normalize_text(value)
        if value and value not in output:
            output.append(value)
    return output


def explicit_species_candidates(text: str) -> list[dict[str, Any]]:
    """Find spans in structurally organism-like positions such as in/from/for."""
    words = normalize_text(text).split()
    cue_words = {"in", "from", "for", "among", "within", "of", "with"}
    candidates: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()

    def add(phrase: str, cue: str, distance: int, width: int) -> None:
        phrase = normalize_text(phrase)
        key = (phrase, cue)
        if not phrase or key in seen:
            return
        seen.add(key)
        candidates.append({"phrase": phrase, "cue": cue, "distance": distance, "width": width})

    for index, token in enumerate(words):
        if token not in cue_words:
            continue
        tail = words[index + 1 : min(len(words), index + 6)]
        for distance in range(len(tail)):
            add(tail[distance], token, distance, 1)
            if distance + 1 < len(tail):
                add(" ".join(tail[distance : distance + 2]), token, distance, 2)

    candidates.sort(key=lambda row: (row["distance"], -row["width"], row["phrase"]))
    return candidates


def compact_query_species_candidates(text: str) -> list[dict[str, str]]:
    """Candidates for compact forms like 'human TP53'; later protein-validated."""
    words = normalize_text(text).split()
    if len(words) < 2 or len(words) > 4:
        return []

    output: list[dict[str, str]] = []
    seen: set[str] = set()

    def add(value: str, side: str) -> None:
        value = normalize_text(value)
        if value and value not in seen:
            seen.add(value)
            output.append({"phrase": value, "side": side})

    add(words[0], "left")
    add(words[-1], "right")
    if len(words) >= 3:
        add(" ".join(words[:2]), "left")
        add(" ".join(words[-2:]), "right")
    return output


def candidate_species_phrases(text: str) -> list[str]:
    """Backward-compatible helper used by tests and diagnostics."""
    values: list[str] = []
    for candidate in [*explicit_species_candidates(text), *compact_query_species_candidates(text)]:
        phrase = candidate["phrase"]
        for value in _lookup_variants(phrase):
            if value not in values:
                values.append(value)
    return values[:16]


def _record_exact_match(record: dict[str, Any], lookup_term: str) -> tuple[int, int, str] | None:
    if taxonomy_record_is_placeholder(record):
        return None
    target = normalize_text(lookup_term)
    rank = normalize_text(record.get("rank", ""))
    rank_score = 0 if rank in SPECIFIC_RANKS else 1

    best: tuple[int, int, str] | None = None
    for kind, name in taxonomy_names(record):
        candidate = normalize_text(name)
        if candidate != target:
            continue
        kind_score = {"scientific name": 0, "common name": 1, "synonym": 2}.get(kind, 3)
        score = (rank_score, kind_score, candidate)
        if best is None or score < best:
            best = score
    return best


def _matches_for_candidate(
    candidate: dict[str, Any],
    rows_by_term: dict[str, list[dict[str, Any]]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    specific: dict[Any, dict[str, Any]] = {}
    broad: dict[Any, dict[str, Any]] = {}

    for lookup_priority, lookup_term in enumerate(_lookup_variants(candidate["phrase"])):
        for record in rows_by_term.get(lookup_term, []):
            score = _record_exact_match(record, lookup_term)
            if score is None:
                continue
            item = {
                "input_phrase": candidate["phrase"],  # preserve what was actually in the sentence
                "lookup_term": lookup_term,
                "record": record,
                "priority": (candidate.get("distance", 0), lookup_priority, *score),
            }
            taxon_id = record.get("taxonId")
            bucket = specific if normalize_text(record.get("rank", "")) in SPECIFIC_RANKS else broad
            existing = bucket.get(taxon_id)
            if existing is None or item["priority"] < existing["priority"]:
                bucket[taxon_id] = item

    return (
        sorted(specific.values(), key=lambda item: item["priority"]),
        sorted(broad.values(), key=lambda item: item["priority"]),
    )


def _remove_candidate(text: str, phrase: str) -> str:
    value = normalize_text(text)
    target = normalize_text(phrase)
    value = re.sub(rf"(?<![a-z0-9]){re.escape(target)}(?![a-z0-9])", " ", value, count=1)
    value = re.sub(r"\b(?:in|from|for|of|the|a|an)\b", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def _compact_candidate_has_direct_protein(
    text: str,
    phrase: str,
    taxonomy_id: int | str,
    client: DataClient,
) -> bool:
    remainder = _remove_candidate(text, phrase)
    if not remainder:
        return False
    try:
        rows = client.uniprot_search(f"({remainder}) AND (organism_id:{taxonomy_id})", size=30)
    except SourceUnavailable:
        return False
    return any(direct_match_details(record, remainder) for record in rows)


def resolve_species(text: str, client: DataClient) -> dict[str, Any]:
    """Resolve only structurally supported organisms; taxonomy validates, never invents."""
    explicit = explicit_species_candidates(text)
    warnings: list[dict[str, str]] = []
    broad_fallback: list[dict[str, Any]] = []

    if explicit:
        terms: list[str] = []
        for candidate in explicit:
            for value in _lookup_variants(candidate["phrase"]):
                if value not in terms:
                    terms.append(value)
        rows_by_term, search_warnings = client.concurrent_taxonomy_search(terms)
        warnings.extend(search_warnings)

        # Work through candidates in grammar priority order. Literal/singular
        # lookup terms may differ, but input_phrase remains the original wording
        # (e.g. "humans"), so the query parser can remove it correctly.
        for candidate in explicit:
            specific, broad = _matches_for_candidate(candidate, rows_by_term)
            if len(specific) == 1:
                return {"status": "resolved", "match": specific[0], "warnings": warnings}
            if len(specific) > 1:
                strongest = specific[0]["priority"][:-1]
                tied = [item for item in specific if item["priority"][:-1] == strongest]
                if len(tied) == 1:
                    return {"status": "resolved", "match": tied[0], "warnings": warnings}
                return {"status": "ambiguous", "matches": tied, "warnings": warnings}
            if broad and not broad_fallback:
                broad_fallback = broad[:8]

        if broad_fallback:
            return {"status": "broad_only", "matches": broad_fallback, "warnings": warnings}

        return {"status": "not_found", "matches": [], "warnings": warnings}

    # Compact syntax such as "human TP53" gets an organism interpretation
    # only when the remainder is a direct protein/gene identity match there.
    for candidate in compact_query_species_candidates(text):
        terms = _lookup_variants(candidate["phrase"])
        rows_by_term, search_warnings = client.concurrent_taxonomy_search(terms)
        warnings.extend(search_warnings)
        specific, _ = _matches_for_candidate({**candidate, "distance": 0}, rows_by_term)
        validated = [
            item
            for item in specific
            if item["record"].get("taxonId") is not None
            and _compact_candidate_has_direct_protein(
                text,
                candidate["phrase"],
                item["record"].get("taxonId"),
                client,
            )
        ]
        if len(validated) == 1:
            return {"status": "resolved", "match": validated[0], "warnings": warnings}
        if len(validated) > 1:
            return {"status": "ambiguous", "matches": validated, "warnings": warnings}

    return {"status": "not_found", "matches": [], "warnings": warnings}
