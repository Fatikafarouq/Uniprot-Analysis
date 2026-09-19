from __future__ import annotations

from collections import defaultdict
from typing import Any

from .evidence import annotation_match_contexts
from .http import DataClient, SourceUnavailable
from .records import direct_match_details, get_gene_name, get_protein_name
from .text import normalize_text


def _quote_query(value: str) -> str:
    return '"' + value.replace('"', "") + '"'


def _scope(query: str, taxon_id: int | str | None) -> str:
    if taxon_id:
        return f"({query}) AND (organism_id:{taxon_id})"
    return f"({query})"


def direct_search(query: str, client: DataClient, taxon_id: int | str | None = None) -> dict[str, Any]:
    """Run independent search variants concurrently, then keep only literal identity matches."""
    normalized = normalize_text(query)
    variants = [
        _scope(f"gene_exact:{_quote_query(normalized)}", taxon_id),
        _scope(f"gene:{_quote_query(normalized)}", taxon_id),
        _scope(_quote_query(normalized), taxon_id),
    ]
    rows_by_query, warnings = client.concurrent_uniprot_search(variants, size=100)

    dedup: dict[str, dict[str, Any]] = {}
    for rows in rows_by_query.values():
        for record in rows:
            accession = record.get("primaryAccession")
            if accession:
                dedup[str(accession)] = record

    matches: list[dict[str, Any]] = []
    for record in dedup.values():
        details = direct_match_details(record, normalized)
        if details:
            matches.append({"record": record, "match": details})

    return {"matches": matches, "warnings": warnings}


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


def fetch_all_records_for_gene(gene: str, taxonomy_id: int | str, client: DataClient) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
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


def discovery_search(
    query: str,
    client: DataClient,
    taxon_id: int | str | None = None,
    limit: int = 20,
) -> dict[str, Any]:
    """Second-screen discovery. Results require traceable evidence in one UniProt-exposed field."""
    search_query = _scope(_quote_query(normalize_text(query)), taxon_id)
    try:
        rows = client.uniprot_search(search_query, size=100)
        warnings: list[dict[str, str]] = []
    except SourceUnavailable as exc:
        rows = []
        warnings = [exc.as_warning()]

    output: list[dict[str, Any]] = []
    for record in rows:
        contexts = annotation_match_contexts(record, query)
        if not contexts:
            continue
        organism = record.get("organism", {}) or {}
        output.append(
            {
                "accession": record.get("primaryAccession"),
                "gene": get_gene_name(record),
                "protein_name": get_protein_name(record),
                "organism": organism.get("commonName") or organism.get("scientificName"),
                "taxon_id": organism.get("taxonId"),
                "why": contexts,
            }
        )
        if len(output) >= limit:
            break
    return {"results": output, "warnings": warnings}
