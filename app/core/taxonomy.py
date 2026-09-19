from __future__ import annotations

from collections import defaultdict
from typing import Any

from .http import DataClient
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


def candidate_species_phrases(text: str) -> list[str]:
    """Small, conservative candidate set. No WordNet and no species alias table."""
    normalized = normalize_text(text)
    if not normalized:
        return []

    tokens = normalized.split()
    candidates: list[str] = []

    def add(value: str) -> None:
        value = normalize_text(value)
        if value and value not in candidates:
            candidates.append(value)

    # Strongest signal: wording after common organism-introducing prepositions.
    padded = f" {normalized} "
    for marker in (" in ", " from ", " of "):
        if marker in padded:
            tail = padded.rsplit(marker, 1)[-1].strip()
            tail_tokens = tail.split()
            for width in range(min(4, len(tail_tokens)), 0, -1):
                add(" ".join(tail_tokens[-width:]))

    # Also consider short suffixes only; this avoids firing taxonomy lookup for every term.
    for width in range(min(3, len(tokens)), 0, -1):
        add(" ".join(tokens[-width:]))

    for value in list(candidates):
        singular = singularize_phrase(value)
        if singular != value:
            add(singular)

    return candidates[:8]


def _match_priority(phrase: str, record: dict[str, Any]) -> tuple[int, int, str] | None:
    target = normalize_text(phrase)
    if not target or taxonomy_record_is_placeholder(record):
        return None

    best: tuple[int, int, str] | None = None
    for kind, name in taxonomy_names(record):
        candidate = normalize_text(name)
        if not candidate:
            continue
        if candidate == target:
            kind_score = {"scientific name": 0, "common name": 1, "synonym": 2}.get(kind, 3)
            rank_score = 0 if str(record.get("rank", "")).lower() in SPECIFIC_RANKS else 1
            score = (rank_score, kind_score, candidate)
            if best is None or score < best:
                best = score
    return best


def resolve_species(text: str, client: DataClient) -> dict[str, Any]:
    phrases = candidate_species_phrases(text)
    if not phrases:
        return {"status": "not_found", "matches": [], "warnings": []}

    rows_by_phrase, warnings = client.concurrent_taxonomy_search(phrases)
    matches: list[dict[str, Any]] = []

    for phrase in phrases:
        for record in rows_by_phrase.get(phrase, []):
            priority = _match_priority(phrase, record)
            if priority is None:
                continue
            matches.append(
                {
                    "input_phrase": phrase,
                    "record": record,
                    "priority": priority,
                }
            )

    # Deduplicate by taxon ID and preserve the strongest literal match.
    best_by_taxon: dict[Any, dict[str, Any]] = {}
    for match in matches:
        taxon_id = match["record"].get("taxonId")
        existing = best_by_taxon.get(taxon_id)
        if existing is None or match["priority"] < existing["priority"]:
            best_by_taxon[taxon_id] = match

    specific = [
        match
        for match in best_by_taxon.values()
        if str(match["record"].get("rank", "")).lower() in SPECIFIC_RANKS
    ]
    specific.sort(key=lambda item: item["priority"])

    if len(specific) == 1:
        return {"status": "resolved", "match": specific[0], "warnings": warnings}
    if len(specific) > 1:
        strongest = specific[0]["priority"][:2]
        tied = [item for item in specific if item["priority"][:2] == strongest]
        if len(tied) == 1:
            return {"status": "resolved", "match": tied[0], "warnings": warnings}
        return {"status": "ambiguous", "matches": tied, "warnings": warnings}

    broad = list(best_by_taxon.values())
    if broad:
        broad.sort(key=lambda item: item["priority"])
        return {"status": "broad_only", "matches": broad[:8], "warnings": warnings}

    return {"status": "not_found", "matches": [], "warnings": warnings}
