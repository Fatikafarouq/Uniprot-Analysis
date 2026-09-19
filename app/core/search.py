from __future__ import annotations

from typing import Any, Iterable

from .evidence import annotation_match_contexts_for_phrases, text_matches_query
from .http import DataClient, SourceUnavailable
from .records import direct_match_details, get_gene_name, get_protein_name
from .text import normalize_text


def _quote_query(value: str) -> str:
    return '"' + value.replace('"', '') + '"'


def _scope(query: str, taxon_id: int | str | None) -> str:
    if taxon_id:
        return f"({query}) AND (organism_id:{taxon_id})"
    return f"({query})"


def _dedupe_records(rows_by_query: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    dedup: dict[str, dict[str, Any]] = {}
    anonymous = 0
    for rows in rows_by_query.values():
        for record in rows:
            accession = record.get("primaryAccession")
            if accession:
                dedup[str(accession)] = record
            else:
                anonymous += 1
                dedup[f"__anonymous_{anonymous}"] = record
    return list(dedup.values())


def direct_search(query: str, client: DataClient, taxon_id: int | str | None = None) -> dict[str, Any]:
    """Run independent search variants concurrently, then keep only literal identity matches."""
    normalized = normalize_text(query)
    variants = [
        _scope(f"gene_exact:{_quote_query(normalized)}", taxon_id),
        _scope(f"gene:{_quote_query(normalized)}", taxon_id),
        _scope(_quote_query(normalized), taxon_id),
    ]
    rows_by_query, warnings = client.concurrent_uniprot_search(variants, size=100)

    matches: list[dict[str, Any]] = []
    for record in _dedupe_records(rows_by_query):
        details = direct_match_details(record, normalized)
        if details:
            matches.append({"record": record, "match": details, "matched_phrase": normalized})

    return {"matches": matches, "warnings": warnings}


def direct_search_phrases(
    phrases: Iterable[str], client: DataClient, taxon_id: int | str | None = None
) -> dict[str, Any]:
    merged: dict[str, dict[str, Any]] = {}
    warnings: list[dict[str, str]] = []

    for phrase in dict.fromkeys(normalize_text(value) for value in phrases if normalize_text(value)):
        payload = direct_search(phrase, client, taxon_id=taxon_id)
        warnings.extend(payload.get("warnings", []))
        for item in payload.get("matches", []):
            accession = item["record"].get("primaryAccession") or f"anon-{id(item['record'])}"
            merged[str(accession)] = item

    return {"matches": list(merged.values()), "warnings": _dedupe_warnings(warnings)}


def group_direct_matches(matches: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, Any], dict[str, Any]] = {}
    for item in matches:
        record = item["record"]
        organism = record.get("organism", {}) or {}
        gene = get_gene_name(record) or item["match"].get("gene") or get_protein_name(record) or "Unknown"
        taxon_id = organism.get("taxonId")
        key = (str(gene), taxon_id)
        group = groups.setdefault(
            key,
            {
                "gene": gene,
                "taxon_id": taxon_id,
                "organism": organism.get("commonName") or organism.get("scientificName") or "Unknown organism",
                "scientific_name": organism.get("scientificName"),
                "record_count_in_search": 0,
                "accessions": [],
            },
        )
        group["record_count_in_search"] += 1
        accession = record.get("primaryAccession")
        if accession and accession not in group["accessions"]:
            group["accessions"].append(accession)
    return sorted(groups.values(), key=lambda row: (str(row["organism"]), str(row["gene"])))


def fetch_all_records_for_gene(
    gene: str, taxonomy_id: int | str, client: DataClient
) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """One streamed UniProt request, not one request per accession."""
    warnings: list[dict[str, str]] = []
    exact = f"(gene_exact:{gene}) AND (organism_id:{taxonomy_id})"
    try:
        records = client.uniprot_stream(exact)
    except SourceUnavailable as exc:
        warnings.append(exc.as_warning())
        records = []

    if not records:
        fallback = f"(gene:{gene}) AND (organism_id:{taxonomy_id})"
        try:
            records = client.uniprot_stream(fallback)
        except SourceUnavailable as exc:
            warnings.append(exc.as_warning())
            records = []

    return records, warnings


def _free_text_query(phrase: str, scope_clause: str | None = None) -> str:
    # Mirrors the successful Colab behavior: preserve the phrase and let
    # UniProt search its indexed text. Evidence filtering happens afterwards.
    if scope_clause:
        return f"({phrase}) AND ({scope_clause})"
    return f"({phrase})"


def search_uniprot_variants(
    phrases: Iterable[str],
    client: DataClient,
    scope_clause: str | None = None,
    size: int = 100,
) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    clean = list(dict.fromkeys(normalize_text(value) for value in phrases if normalize_text(value)))
    queries = [_free_text_query(phrase, scope_clause) for phrase in clean]
    if not queries:
        return [], []
    rows_by_query, warnings = client.concurrent_uniprot_search(queries, size=size)
    return _dedupe_records(rows_by_query), warnings


def _record_source(record: dict[str, Any]) -> dict[str, Any]:
    organism = record.get("organism", {}) or {}
    scientific = organism.get("scientificName")
    common = organism.get("commonName")
    return {
        "friendly": common or scientific or "Unknown organism",
        "scientific": scientific,
        "common": common,
        "taxon_id": organism.get("taxonId"),
    }


def _contexts_with_source_name(
    record: dict[str, Any], phrases: list[str], max_items: int = 3
) -> tuple[list[dict[str, str]], str | None]:
    source = _record_source(record)
    for phrase in phrases:
        contexts: list[dict[str, str]] = []
        for label, value in (
            ("Source organism scientific name", source.get("scientific")),
            ("Source organism common name", source.get("common")),
        ):
            if value and text_matches_query(str(value), phrase):
                contexts.append({"source": label, "text": str(value)})

        annotation_contexts, _ = annotation_match_contexts_for_phrases(record, [phrase], max_items=max_items)
        for context in annotation_contexts:
            signature = (context.get("source"), context.get("text"))
            if any((item.get("source"), item.get("text")) == signature for item in contexts):
                continue
            contexts.append(context)
            if len(contexts) >= max_items:
                break

        if contexts:
            return contexts[:max_items], phrase
    return [], None


def _discovery_item(
    record: dict[str, Any],
    contexts: list[dict[str, str]],
    matched_phrase: str | None,
    *,
    bucket: str,
    relationship: str,
    host: dict[str, Any] | None = None,
) -> dict[str, Any]:
    organism = record.get("organism", {}) or {}
    return {
        "accession": record.get("primaryAccession"),
        "gene": get_gene_name(record),
        "protein_name": get_protein_name(record),
        "organism": organism.get("commonName") or organism.get("scientificName"),
        "scientific_name": organism.get("scientificName"),
        "taxon_id": organism.get("taxonId"),
        "why": contexts,
        "matched_phrase": matched_phrase,
        "bucket": bucket,
        "relationship": relationship,
        "host": host,
    }


def _matching_host(record: dict[str, Any], taxonomy_id: int | str) -> dict[str, Any] | None:
    for host in record.get("organismHosts", []) or []:
        if str(host.get("taxonId")) == str(taxonomy_id):
            return host
    return None


def _dedupe_warnings(warnings: list[dict[str, str]]) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for warning in warnings:
        key = (str(warning.get("source", "")), str(warning.get("message", "")))
        if key in seen:
            continue
        seen.add(key)
        output.append(warning)
    return output


def discovery_search(
    query: str | Iterable[str],
    client: DataClient,
    taxon_id: int | str | None = None,
    organism_name: str | None = None,
    limit: int = 24,
) -> dict[str, Any]:
    """Colab-style discovery: mentioned organism, virus-host, then global evidence."""
    if isinstance(query, str):
        phrases = [normalize_text(query)]
    else:
        phrases = [normalize_text(value) for value in query]
    phrases = list(dict.fromkeys(value for value in phrases if value))

    warnings: list[dict[str, str]] = []
    output: list[dict[str, Any]] = []
    seen: set[str] = set()

    def append_item(item: dict[str, Any]) -> None:
        accession = str(item.get("accession") or "")
        if not accession or accession in seen or len(output) >= limit:
            return
        seen.add(accession)
        output.append(item)

    if taxon_id is not None:
        source_rows, source_warnings = search_uniprot_variants(
            phrases, client, scope_clause=f"organism_id:{taxon_id}", size=100
        )
        warnings.extend(source_warnings)
        for record in source_rows:
            contexts, matched_phrase = annotation_match_contexts_for_phrases(record, phrases, max_items=3)
            if contexts:
                append_item(
                    _discovery_item(
                        record,
                        contexts,
                        matched_phrase,
                        bucket="mentioned_organism",
                        relationship=f"Protein from {organism_name or 'the organism you mentioned'}",
                    )
                )

        virus_rows, virus_warnings = search_uniprot_variants(
            phrases, client, scope_clause=f"virus_host_id:{taxon_id}", size=100
        )
        warnings.extend(virus_warnings)
        for record in virus_rows:
            contexts, matched_phrase = _contexts_with_source_name(record, phrases, max_items=3)
            if not contexts:
                continue
            append_item(
                _discovery_item(
                    record,
                    contexts,
                    matched_phrase,
                    bucket="virus_host",
                    relationship=(
                        f"Viral protein from a virus UniProt lists with "
                        f"{organism_name or 'the mentioned organism'} as a host"
                    ),
                    host=_matching_host(record, taxon_id),
                )
            )

    global_rows, global_warnings = search_uniprot_variants(phrases, client, scope_clause=None, size=100)
    warnings.extend(global_warnings)
    for record in global_rows:
        source = _record_source(record)
        # Mentioned-organism records already have a clearer bucket above.
        if taxon_id is not None and str(source.get("taxon_id")) == str(taxon_id):
            continue
        contexts, matched_phrase = _contexts_with_source_name(record, phrases, max_items=3)
        if not contexts:
            continue
        append_item(
            _discovery_item(
                record,
                contexts,
                matched_phrase,
                bucket="global",
                relationship=(
                    "Global UniProt discovery"
                    if taxon_id is None
                    else "Protein from another UniProt source organism"
                ),
            )
        )

    return {
        "results": output,
        "warnings": _dedupe_warnings(warnings),
        "phrases": phrases,
    }
