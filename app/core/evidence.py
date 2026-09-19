from __future__ import annotations

from typing import Any, Iterable

from .records import get_all_gene_labels, get_all_protein_search_names
from .text import meaningful_tokens, normalize_text

SEMANTIC_XREF_KEY_WORDS = {"term", "pathway", "disease", "phenotype", "process", "function"}


def literature_title_identifies_record(record: dict[str, Any], title: str) -> bool:
    title_norm = normalize_text(title)
    names = [*get_all_gene_labels(record), *get_all_protein_search_names(record)]
    return any(normalize_text(name) in title_norm for name in names if normalize_text(name))


def iter_discovery_strings(record: dict[str, Any]):
    for comment in record.get("comments", []) or []:
        comment_type = str(comment.get("commentType") or "Comment")
        for text in comment.get("texts", []) or []:
            if isinstance(text, dict) and text.get("value"):
                yield f"UniProt {comment_type} comment", str(text["value"])

    for keyword in record.get("keywords", []) or []:
        if isinstance(keyword, dict) and keyword.get("name"):
            yield "UniProt keyword", str(keyword["name"])

    for feature in record.get("features", []) or []:
        description = feature.get("description")
        if description:
            yield f"UniProt {feature.get('type', 'feature')} annotation", str(description)

    for xref in record.get("uniProtKBCrossReferences", []) or []:
        database = str(xref.get("database") or "database")
        for prop in xref.get("properties", []) or []:
            if not isinstance(prop, dict):
                continue
            key = str(prop.get("key") or "")
            value = prop.get("value")
            normalized_key = normalize_text(key).replace(" ", "")
            if value and any(word in normalized_key for word in SEMANTIC_XREF_KEY_WORDS):
                yield f"{database} cross-reference ({key})", str(value)

    for reference in record.get("references", []) or []:
        citation = reference.get("citation", {}) or {}
        title = citation.get("title")
        if title and literature_title_identifies_record(record, str(title)):
            yield "Reference title linked by UniProt", str(title)


def text_matches_query(text: str, query: str) -> bool:
    haystack = normalize_text(text)
    query_norm = normalize_text(query)
    if not haystack or not query_norm:
        return False
    if query_norm in haystack:
        return True
    tokens = meaningful_tokens(query_norm)
    return bool(tokens) and all(token in haystack for token in tokens)


def annotation_match_contexts(record: dict[str, Any], query: str, max_items: int = 3) -> list[dict[str, str]]:
    """A single source field must support the complete query; fields are never combined."""
    contexts: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for source, value in iter_discovery_strings(record):
        clean = " ".join(str(value).split())
        if not clean or not text_matches_query(clean, query):
            continue
        signature = (source, clean)
        if signature in seen:
            continue
        seen.add(signature)
        contexts.append({"source": source, "text": clean})
        if len(contexts) >= max_items:
            break
    return contexts


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

    return {
        "summary": summary,
        "category_counts": category_counts,
        "statements": statements,
    }
