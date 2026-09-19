from __future__ import annotations

import difflib
from typing import Any, Iterable

from .evidence import (
    DISCOVERY_EVIDENCE_TYPES,
    annotation_match_contexts,
    annotation_match_contexts_for_phrases,
    classify_discovery_evidence,
    beginner_connection_explanation,
    recursive_strings,
    shorten_match_text,
    text_matches_query,
)
from .http import DataClient, SourceUnavailable
from .records import (
    direct_match_details,
    get_alternative_names,
    get_function_notes,
    get_gene_name,
    get_protein_name,
)
from .text import normalize_text


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


def _free_text_query(phrase: str, scope_clause: str | None = None) -> str:
    if scope_clause:
        return f"({phrase}) AND ({scope_clause})"
    return f"({phrase})"


def search_uniprot_variants(
    phrases: Iterable[str],
    client: DataClient,
    scope_clause: str | None = None,
    size: int = 100,
) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """Colab-v14 free-text search, but variants are executed concurrently."""
    clean = list(dict.fromkeys(normalize_text(value) for value in phrases if normalize_text(value)))
    queries = [_free_text_query(phrase, scope_clause) for phrase in clean]
    if not queries:
        return [], []
    rows_by_query, warnings = client.concurrent_uniprot_search(queries, size=size)
    return _dedupe_records(rows_by_query), _dedupe_warnings(warnings)


def build_direct_gene_options(initial_results: list[dict[str, Any]], phrase: str) -> dict[str, dict[str, Any]]:
    """Port of the Colab v14 direct-name grouping rule."""
    options: dict[str, dict[str, Any]] = {}
    for record in initial_results:
        gene = get_gene_name(record)
        if not gene:
            continue
        details = direct_match_details(record, phrase)
        if not details:
            continue
        option = options.setdefault(gene, {"search_records": [], "reasons": []})
        option["search_records"].append(record)
        option["reasons"].extend(details.get("reasons", []))
    for info in options.values():
        info["reasons"] = list(dict.fromkeys(info["reasons"]))
    return options


def direct_search_phrases(
    phrases: Iterable[str], client: DataClient, taxon_id: int | str | None = None
) -> dict[str, Any]:
    """Use the same free-text retrieval + name-field filtering as the Colab prototype."""
    clean = list(dict.fromkeys(normalize_text(value) for value in phrases if normalize_text(value)))
    scope = f"organism_id:{taxon_id}" if taxon_id is not None else None
    rows, warnings = search_uniprot_variants(clean, client, scope_clause=scope, size=100)

    matches: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for phrase in clean:
        for record in rows:
            details = direct_match_details(record, phrase)
            if not details:
                continue
            accession = str(record.get("primaryAccession") or id(record))
            key = (accession, phrase)
            if key in seen:
                continue
            seen.add(key)
            matches.append({"record": record, "match": details, "matched_phrase": phrase})

    # Identity-first rule: if any record exactly matches a UniProt gene label
    # or protein name, do not present weaker substring/name-containing hits as
    # competing direct identities. Those weaker hits remain discoverable when
    # there is no exact identity, preserving the Colab behaviour for searches
    # such as "anthrax" -> Anthrax toxin receptor 1/2.
    if any(item.get("match", {}).get("match_strength") == "exact" for item in matches):
        matches = [
            item for item in matches
            if item.get("match", {}).get("match_strength") == "exact"
        ]

    return {"matches": matches, "warnings": warnings, "search_records": rows}


def group_direct_matches(matches: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, Any], dict[str, Any]] = {}
    for item in matches:
        record = item["record"]
        organism = record.get("organism", {}) or {}
        details = item.get("match", {}) or {}
        gene = get_gene_name(record) or details.get("gene") or get_protein_name(record) or "Unknown"
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
                "reasons": [],
            },
        )
        group["record_count_in_search"] += 1
        accession = record.get("primaryAccession")
        if accession and accession not in group["accessions"]:
            group["accessions"].append(accession)
        group["reasons"].extend(details.get("reasons", []))
    for group in groups.values():
        group["reasons"] = list(dict.fromkeys(group["reasons"]))
    return sorted(groups.values(), key=lambda row: (str(row["organism"]), str(row["gene"])))


def fetch_all_records_for_gene(
    gene: str, taxonomy_id: int | str, client: DataClient
) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """Retrieve the whole gene/species set via UniProt's stream endpoint."""
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
    return records, _dedupe_warnings(warnings)


def record_source_organism(record: dict[str, Any]) -> dict[str, Any]:
    organism = record.get("organism", {}) or {}
    scientific = organism.get("scientificName")
    common = organism.get("commonName")
    return {
        "friendly": common or scientific or "Unknown organism",
        "scientific": scientific,
        "common": common,
        "taxon_id": organism.get("taxonId"),
    }


def summarise_option_records(records: list[dict[str, Any]], fallback_gene: str) -> dict[str, Any]:
    names: list[str] = []
    aliases: list[str] = []
    functions: list[str] = []
    for record in records:
        name = get_protein_name(record)
        if name:
            names.append(name)
        aliases.extend(get_alternative_names(record))
        functions.extend(get_function_notes(record))
    return {
        "name": list(dict.fromkeys(names))[0] if names else fallback_gene,
        "aliases": list(dict.fromkeys(aliases)),
        "functions": list(dict.fromkeys(functions)),
    }


def brief_function_from_records(records: list[dict[str, Any]], matched_phrase: str | None, limit: int = 240) -> str | None:
    seen: set[str] = set()
    for record in records:
        accession = str(record.get("primaryAccession") or id(record))
        if accession in seen:
            continue
        seen.add(accession)
        notes = get_function_notes(record)
        if notes:
            return shorten_match_text(notes[0], matched_phrase or "", limit=limit)
    return None


def build_discovery_gene_options(
    initial_results: list[dict[str, Any]], phrases: list[str] | str, limit: int = 12
) -> dict[str, dict[str, Any]]:
    """Colab v14 behavior: group mentioned-organism discovery by gene, not accession."""
    if isinstance(phrases, str):
        phrases = [phrases]
    options: dict[str, dict[str, Any]] = {}
    for record in initial_results:
        gene = get_gene_name(record)
        if not gene:
            continue
        contexts, matched_phrase = annotation_match_contexts_for_phrases(record, phrases, max_items=3)
        if not contexts:
            continue
        if gene not in options:
            if len(options) >= limit:
                continue
            options[gene] = {"search_records": [], "contexts": [], "matched_phrase": matched_phrase}
        option = options[gene]
        option["search_records"].append(record)
        option["contexts"].extend(contexts)
        if option.get("matched_phrase") is None:
            option["matched_phrase"] = matched_phrase

    for info in options.values():
        deduped: list[dict[str, str]] = []
        seen: set[tuple[str | None, str | None]] = set()
        for context in info["contexts"]:
            key = (context.get("source"), context.get("text"))
            if key in seen:
                continue
            seen.add(key)
            deduped.append(context)
        info["contexts"] = deduped[:3]
    return options


def _contexts_with_source_name(
    record: dict[str, Any], phrases: list[str], max_items: int = 3
) -> tuple[list[dict[str, str]], str | None]:
    source = record_source_organism(record)
    for phrase in phrases:
        contexts: list[dict[str, str]] = []
        for label, value in (
            ("Source organism scientific name", source.get("scientific")),
            ("Source organism common name", source.get("common")),
        ):
            if value and text_matches_query(str(value), phrase):
                contexts.append({"source": label, "text": str(value)})

        for context in annotation_match_contexts(record, phrase, max_items=max_items):
            signature = (context.get("source"), context.get("text"))
            if any((item.get("source"), item.get("text")) == signature for item in contexts):
                continue
            contexts.append(context)
            if len(contexts) >= max_items:
                break
        if contexts:
            return contexts[:max_items], phrase
    return [], None


def _matching_host(record: dict[str, Any], taxonomy_id: int | str) -> dict[str, Any] | None:
    for host in record.get("organismHosts", []) or []:
        if str(host.get("taxonId")) == str(taxonomy_id):
            return host
    return None


def build_virus_host_options(
    records: list[dict[str, Any]],
    phrases: list[str],
    taxonomy_id: int | str,
    host_friendly: str | None,
    limit: int = 12,
) -> list[dict[str, Any]]:
    options: list[dict[str, Any]] = []
    seen_accessions: set[str] = set()
    for record in records:
        accession = record.get("primaryAccession")
        if not accession or str(accession) in seen_accessions:
            continue
        contexts, matched_phrase = _contexts_with_source_name(record, phrases, max_items=3)
        if not contexts:
            continue
        source = record_source_organism(record)
        options.append(
            {
                "kind": "virus_host",
                "accession": accession,
                "gene": get_gene_name(record),
                "protein_name": get_protein_name(record) or accession,
                "organism": source["friendly"],
                "scientific_name": source.get("scientific"),
                "taxon_id": source.get("taxon_id"),
                "why": contexts,
                "matched_phrase": matched_phrase,
                "bucket": "virus_host",
                "relationship": f"Viral protein from a virus UniProt lists with {host_friendly or 'the mentioned organism'} as a host",
                "host": _matching_host(record, taxonomy_id),
                "evidence_type": classify_discovery_evidence(contexts),
                "connection_explanation": beginner_connection_explanation(
                    classify_discovery_evidence(contexts), matched_phrase
                ),
                "function_note": brief_function_from_records([record], matched_phrase),
                "inspect": {
                    "mode": "gene" if get_gene_name(record) and source.get("taxon_id") else "accession",
                    "gene": get_gene_name(record),
                    "taxon_id": source.get("taxon_id"),
                    "species_name": source["friendly"],
                    "accession": accession,
                },
            }
        )
        seen_accessions.add(str(accession))
        if len(options) >= limit:
            break
    return options


def build_global_discovery_options(
    records: list[dict[str, Any]],
    phrases: list[str],
    mentioned_taxonomy_id: int | str | None = None,
    excluded_accessions: set[str] | None = None,
    limit: int = 12,
) -> list[dict[str, Any]]:
    excluded_accessions = set(excluded_accessions or [])
    options: list[dict[str, Any]] = []
    seen = set(excluded_accessions)
    for record in records:
        accession = record.get("primaryAccession")
        if not accession or str(accession) in seen:
            continue
        source = record_source_organism(record)
        if mentioned_taxonomy_id is not None and str(source.get("taxon_id")) == str(mentioned_taxonomy_id):
            continue
        contexts, matched_phrase = _contexts_with_source_name(record, phrases, max_items=3)
        if not contexts:
            continue
        gene = get_gene_name(record)
        options.append(
            {
                "kind": "global_discovery",
                "accession": accession,
                "gene": gene,
                "protein_name": get_protein_name(record) or accession,
                "organism": source["friendly"],
                "scientific_name": source.get("scientific"),
                "taxon_id": source.get("taxon_id"),
                "why": contexts,
                "matched_phrase": matched_phrase,
                "bucket": "global",
                "relationship": "Global UniProt discovery" if mentioned_taxonomy_id is None else "Protein from another UniProt source organism",
                "evidence_type": classify_discovery_evidence(contexts),
                "connection_explanation": beginner_connection_explanation(
                    classify_discovery_evidence(contexts), matched_phrase
                ),
                "function_note": brief_function_from_records([record], matched_phrase),
                "inspect": {
                    "mode": "gene" if gene and source.get("taxon_id") else "accession",
                    "gene": gene,
                    "taxon_id": source.get("taxon_id"),
                    "species_name": source["friendly"],
                    "accession": accession,
                },
            }
        )
        seen.add(str(accession))
        if len(options) >= limit:
            break
    return options


def discovery_search(
    query: str | Iterable[str],
    client: DataClient,
    taxon_id: int | str | None = None,
    organism_name: str | None = None,
    limit: int = 36,
    source_rows: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Colab three-path discovery flow with web-oriented request reuse.

    When ``source_rows`` are supplied, they are reused for the mentioned
    organism instead of repeating the same UniProt request that was already
    used for direct-name matching. Virus-host and global branches are fetched
    concurrently so discovery does not wait on them one after another.
    """
    if isinstance(query, str):
        phrases = [normalize_text(query)]
    else:
        phrases = [normalize_text(value) for value in query]
    phrases = list(dict.fromkeys(value for value in phrases if value))

    warnings: list[dict[str, str]] = []
    output: list[dict[str, Any]] = []
    excluded_accessions: set[str] = set()

    if taxon_id is not None:
        if source_rows is None:
            source_rows, source_warnings = search_uniprot_variants(
                phrases, client, scope_clause=f"organism_id:{taxon_id}", size=100
            )
            warnings.extend(source_warnings)

        source_options = build_discovery_gene_options(source_rows or [], phrases, limit=12)
        for gene, info in source_options.items():
            summary = summarise_option_records(info["search_records"], gene)
            first = info["search_records"][0]
            accession = first.get("primaryAccession")
            for row in info["search_records"]:
                if row.get("primaryAccession"):
                    excluded_accessions.add(str(row["primaryAccession"]))
            output.append(
                {
                    "kind": "source_discovery",
                    "accession": accession,
                    "gene": gene,
                    "protein_name": summary["name"],
                    "organism": organism_name or record_source_organism(first)["friendly"],
                    "scientific_name": record_source_organism(first).get("scientific"),
                    "taxon_id": taxon_id,
                    "why": info["contexts"],
                    "matched_phrase": info.get("matched_phrase"),
                    "bucket": "mentioned_organism",
                    "relationship": f"Protein from {organism_name or 'the organism you mentioned'}",
                    "evidence_type": classify_discovery_evidence(info["contexts"]),
                    "connection_explanation": beginner_connection_explanation(
                        classify_discovery_evidence(info["contexts"]), info.get("matched_phrase")
                    ),
                    "evidence_record_count": len(info["search_records"]),
                    "function_note": brief_function_from_records(info["search_records"], info.get("matched_phrase")),
                    "inspect": {
                        "mode": "gene",
                        "gene": gene,
                        "taxon_id": taxon_id,
                        "species_name": organism_name or record_source_organism(first)["friendly"],
                        "accession": accession,
                    },
                }
            )

        # Fetch the two remaining discovery paths in one concurrent batch.
        query_groups: dict[str, list[str]] = {"virus": [], "global": []}
        for phrase in phrases:
            query_groups["virus"].append(_free_text_query(phrase, f"virus_host_id:{taxon_id}"))
            query_groups["global"].append(_free_text_query(phrase, None))

        all_queries = [*query_groups["virus"], *query_groups["global"]]
        rows_by_query, concurrent_warnings = client.concurrent_uniprot_search(all_queries, size=100)
        warnings.extend(concurrent_warnings)

        virus_rows = _dedupe_records({q: rows_by_query.get(q, []) for q in query_groups["virus"]})
        global_rows = _dedupe_records({q: rows_by_query.get(q, []) for q in query_groups["global"]})

        virus_options = build_virus_host_options(
            virus_rows, phrases, taxon_id, organism_name, limit=12
        )
        output.extend(virus_options)
        for item in virus_options:
            if item.get("accession"):
                excluded_accessions.add(str(item["accession"]))
    else:
        global_rows, global_warnings = search_uniprot_variants(
            phrases, client, scope_clause=None, size=100
        )
        warnings.extend(global_warnings)

    global_options = build_global_discovery_options(
        global_rows,
        phrases,
        mentioned_taxonomy_id=taxon_id,
        excluded_accessions=excluded_accessions,
        limit=12,
    )
    output.extend(global_options)

    return {
        "results": output[:limit],
        "warnings": _dedupe_warnings(warnings),
        "phrases": phrases,
        "evidence_types": DISCOVERY_EVIDENCE_TYPES,
    }


def uniprot_spelling_vocabulary(records: list[dict[str, Any]]) -> set[str]:
    vocabulary: set[str] = set()

    def add_text(value: Any) -> None:
        if not value:
            return
        for token in normalize_text(value).split():
            if len(token) >= 4 and any(ch.isalpha() for ch in token):
                vocabulary.add(token)

    for record in records:
        add_text(get_gene_name(record))
        add_text(get_protein_name(record))
        for name in get_alternative_names(record):
            add_text(name)
        for comment in record.get("comments", []) or []:
            for value in recursive_strings(comment):
                add_text(value)
        for keyword in record.get("keywords", []) or []:
            add_text(keyword.get("name") if isinstance(keyword, dict) else keyword)
    return vocabulary


def suggest_uniprot_wording_correction(
    phrase: str, taxonomy_id: int | str, client: DataClient
) -> tuple[dict[str, Any] | None, list[dict[str, str]]]:
    """Database-grounded spelling correction from Colab v14; WordNet is deliberately not used."""
    tokens = normalize_text(phrase).split()
    if not tokens:
        return None, []

    candidates: list[dict[str, Any]] = []
    warnings: list[dict[str, str]] = []

    for index, token in enumerate(tokens):
        if len(token) < 4 or not any(ch.isalpha() for ch in token):
            continue
        anchor_tokens = [value for position, value in enumerate(tokens) if position != index and len(value) >= 3]
        if not anchor_tokens:
            continue
        anchor_phrase = " ".join(anchor_tokens)
        try:
            anchor_results = client.uniprot_search(
                f"({anchor_phrase}) AND (organism_id:{taxonomy_id})", size=100
            )
        except SourceUnavailable as exc:
            warnings.append(exc.as_warning())
            continue
        if not anchor_results:
            continue

        vocabulary = uniprot_spelling_vocabulary(anchor_results)
        close_words = difflib.get_close_matches(token, sorted(vocabulary), n=6, cutoff=0.88)
        for replacement in close_words:
            if replacement == token:
                continue
            ratio = difflib.SequenceMatcher(None, token, replacement).ratio()
            corrected_tokens = list(tokens)
            corrected_tokens[index] = replacement
            corrected_phrase = " ".join(corrected_tokens)
            if corrected_phrase == normalize_text(phrase):
                continue
            try:
                validation_results = client.uniprot_search(
                    f"({corrected_phrase}) AND (organism_id:{taxonomy_id})", size=10
                )
            except SourceUnavailable as exc:
                warnings.append(exc.as_warning())
                continue
            if validation_results:
                candidates.append(
                    {
                        "phrase": corrected_phrase,
                        "ratio": ratio,
                        "changed_from": token,
                        "changed_to": replacement,
                    }
                )

    if not candidates:
        return None, _dedupe_warnings(warnings)
    candidates.sort(key=lambda item: item["ratio"], reverse=True)
    best = candidates[0]
    if len(candidates) > 1:
        second = candidates[1]
        if second["phrase"] != best["phrase"] and abs(second["ratio"] - best["ratio"]) < 0.015:
            return None, _dedupe_warnings(warnings)
    return best, _dedupe_warnings(warnings)
