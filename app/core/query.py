from __future__ import annotations

import re

from .text import normalize_text

# Generic request wording only. These rules are anchored at the LEFT edge so
# biological words inside the concept are never deleted simply to obtain hits.
REQUEST_PREFIX_PATTERNS = [
    r"^(?:please\s+)?(?:give|show|tell|find|list)\s+(?:me\s+)?",
    r"^(?:can|could|would)\s+you\s+(?:please\s+)?(?:give|show|tell|find|list)\s+(?:me\s+)?",
    r"^(?:the\s+)?(?:one|thing|something)\s+(?:that|which)\s+(?:causes?|caused|causing|makes?|made|making|influences?|influenced|influencing|affects?|affected|affecting)\s+",
    r"^(?:what|which)\s+(?:causes?|caused|causing|makes?|made|making|influences?|influenced|influencing|affects?|affected|affecting)\s+",
    r"^(?:what|which)\s+(?:protein|proteins|gene|genes)\s+(?:is|are|does|do|can|could|might|may)?\s*",
    r"^(?:associated|related|linked|connected|involved)\s+(?:with|to|in)\s+",
    r"^(?:most\s+)?important\s+(?:in\s+)?",
    r"^(?:causes?|caused|causing|makes?|made|making|influences?|influenced|influencing|affects?|affected|affecting)\s+",
    r"^susceptible\s+to\s+",
]


def remove_literal_phrase_once(text: str, phrase: str | None) -> str:
    cleaned = normalize_text(text)
    target = normalize_text(phrase)
    if not cleaned or not target:
        return cleaned

    result = re.sub(
        rf"(?<![a-z0-9]){re.escape(target)}(?![a-z0-9])",
        " ",
        cleaned,
        count=1,
    )
    return re.sub(r"\s+", " ", result).strip()


def extract_biological_phrase(text: str, species_phrase: str | None) -> str:
    """Remove only a confirmed organism mention and dangling grammar around it."""
    cleaned = remove_literal_phrase_once(text, species_phrase)

    # When an organism appears at the end ("anthrax in humans"), removing the
    # organism leaves a dangling connector. Remove only edge grammar, never a
    # biological word from the middle of the query.
    cleaned = re.sub(r"\b(?:in|from|for|of|the|a|an)\b\s*$", " ", cleaned)
    cleaned = re.sub(r"^(?:in|from|for|of|the|a|an)\b\s*", " ", cleaned)
    return re.sub(r"\s+", " ", cleaned).strip()


def strip_leading_request_scaffolding(text: str) -> str:
    value = normalize_text(text)
    if not value:
        return value

    for _ in range(8):
        before = value
        for pattern in REQUEST_PREFIX_PATTERNS:
            value = re.sub(pattern, "", value, count=1).strip()
        value = re.sub(r"\s+", " ", value).strip()
        if value == before:
            break
    return value


def concept_orthographic_variants(phrase: str) -> list[str]:
    """Meaning-preserving spacing/hyphen variants; no biological synonym table."""
    phrase = normalize_text(phrase)
    if not phrase:
        return []

    tokens = phrase.split()
    variants: list[str] = []

    def add(value: str) -> None:
        value = normalize_text(value)
        if value and value not in variants:
            variants.append(value)

    add(phrase)
    if 2 <= len(tokens) <= 3:
        add("-".join(tokens))
        add("".join(tokens))
    return variants


def broader_concept_levels(strict_concept: str) -> list[str]:
    """Transparent right-hand fallbacks used only after the full concept fails."""
    tokens = normalize_text(strict_concept).split()
    levels: list[str] = []
    for start in range(1, len(tokens)):
        value = " ".join(tokens[start:]).strip()
        if value and value not in levels:
            levels.append(value)
    return levels


def extract_search_concepts(text: str, species_phrase: str | None) -> dict[str, object]:
    """Parse natural-language wording without guessing biological intent."""
    raw = normalize_text(extract_biological_phrase(text, species_phrase))
    primary = strip_leading_request_scaffolding(raw) or raw
    return {
        "raw": raw,
        "primary": primary,
        "strict_variants": concept_orthographic_variants(primary),
        "related_levels": broader_concept_levels(primary),
    }
