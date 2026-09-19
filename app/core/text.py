from __future__ import annotations

import re


def normalize_text(value: object) -> str:
    text = str(value or "").lower().strip()
    text = re.sub(r"[^a-z0-9\s.\-]", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def strip_version(identifier: str | None) -> str | None:
    if not identifier:
        return None
    return re.sub(r"\.\d+$", "", str(identifier))


def singularize_phrase(phrase: str) -> str:
    words = normalize_text(phrase).split()
    if not words:
        return normalize_text(phrase)

    last = words[-1]
    if last.endswith("ies") and len(last) > 4:
        last = last[:-3] + "y"
    elif last.endswith("s") and not last.endswith("ss") and len(last) > 3:
        last = last[:-1]
    words[-1] = last
    return " ".join(words)


def meaningful_tokens(text: str) -> list[str]:
    stop = {
        "a", "an", "and", "are", "as", "at", "be", "by", "for", "from",
        "in", "is", "of", "on", "or", "protein", "proteins", "gene", "genes",
        "that", "the", "this", "to", "what", "which", "with",
    }
    return [token for token in normalize_text(text).split() if token not in stop and len(token) > 1]
