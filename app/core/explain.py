from __future__ import annotations

import logging
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from typing import Any

from .evidence import function_evidence
from .http import DataClient, SourceUnavailable
from .records import (
    analyse_features,
    existence_kind,
    external_links,
    extract_isoforms,
    get_ensembl_details,
    get_gene_name,
    get_length,
    get_protein_existence,
    get_protein_name,
    get_proteome_ids,
    group_sequence_duplicates,
    review_status,
)
from .text import normalize_text, strip_version

logger = logging.getLogger(__name__)

ISOFORM_EXPLANATION = (
    "A UniProt isoform is an alternative protein sequence described within one UniProtKB entry, "
    "often because of processes such as alternative splicing or alternative promoter usage. "
    "An isoform is not the same thing as a separate UniProt accession."
)


def _gene_centric_accession(data: Any) -> str | None:
    rows = data if isinstance(data, list) else [data] if isinstance(data, dict) else []
    for row in rows:
        for key in ("canonicalProtein", "protein", "accession", "uniprotAccession"):
            value = row.get(key) if isinstance(row, dict) else None
            if isinstance(value, str):
                return value.split("-")[0]
            if isinstance(value, dict):
                nested = value.get("id") or value.get("accession")
                if nested:
                    return str(nested).split("-")[0]
    return None


def _appris_map(data: Any) -> dict[str, str]:
    result: dict[str, str] = {}
    if not isinstance(data, list):
        return result
    for row in data:
        if not isinstance(row, dict):
            continue
        transcript = strip_version(row.get("transcript_id"))
        reliability = row.get("reliability")
        if not transcript or not reliability:
            continue
        existing = result.get(transcript)
        if existing is None or str(reliability).startswith("PRINCIPAL"):
            result[transcript] = str(reliability)
    return result


def get_external_annotation_context(records: list[dict[str, Any]], client: DataClient) -> tuple[dict[str, Any], list[dict[str, str]]]:
    ensembl_gene = None
    proteome_id = None
    for record in records:
        rows = get_ensembl_details(record)
        if not ensembl_gene and rows and rows[0].get("gene"):
            ensembl_gene = strip_version(rows[0]["gene"])
        proteomes = get_proteome_ids(record)
        if not proteome_id and proteomes:
            proteome_id = proteomes[0]

    warnings: list[dict[str, str]] = []
    if not ensembl_gene:
        return {
            "ensembl_gene": None,
            "ensembl_canonical": None,
            "ensembl_species": None,
            "gene_centric_accession": None,
            "appris": {},
        }, warnings

    gene_centric_data = None
    ensembl_data = None

    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = {}
        if proteome_id:
            futures[pool.submit(client.gene_centric, proteome_id, ensembl_gene)] = "gene-centric"
        futures[pool.submit(client.ensembl_lookup, ensembl_gene)] = "Ensembl"

        for future, label in list(futures.items()):
            try:
                value = future.result()
                if label == "gene-centric":
                    gene_centric_data = value
                else:
                    ensembl_data = value
            except SourceUnavailable as exc:
                warnings.append(exc.as_warning())
            except Exception as exc:
                logger.exception("Unexpected external annotation failure: %s", label)
                warnings.append({"source": label, "message": f"{label} was unavailable: {exc}"})

    ensembl_canonical = None
    ensembl_species = None
    if isinstance(ensembl_data, dict):
        ensembl_canonical = ensembl_data.get("canonical_transcript")
        ensembl_species = ensembl_data.get("species")

    appris = {}
    if ensembl_species:
        try:
            appris = _appris_map(client.appris(str(ensembl_species), str(ensembl_gene)))
        except SourceUnavailable as exc:
            warnings.append(exc.as_warning())

    return {
        "ensembl_gene": ensembl_gene,
        "ensembl_canonical": ensembl_canonical,
        "ensembl_species": ensembl_species,
        "gene_centric_accession": _gene_centric_accession(gene_centric_data),
        "appris": appris,
    }, warnings


def build_comparison_records(records: list[dict[str, Any]], external: dict[str, Any]) -> list[dict[str, Any]]:
    duplicates = group_sequence_duplicates(records)
    output: list[dict[str, Any]] = []

    for record in records:
        accession = str(record.get("primaryAccession") or "")
        ensembl_rows = get_ensembl_details(record)
        transcript = ensembl_rows[0].get("transcript") if ensembl_rows else None
        transcript_base = strip_version(transcript)
        canonical = bool(
            transcript_base
            and external.get("ensembl_canonical")
            and transcript_base == strip_version(external.get("ensembl_canonical"))
        )

        output.append(
            {
                "accession": accession,
                "name": get_protein_name(record) or get_gene_name(record) or accession,
                "gene": get_gene_name(record),
                "length": get_length(record),
                "existence": existence_kind(get_protein_existence(record)),
                "review": review_status(record),
                "transcript": transcript,
                "canonical": canonical,
                "appris": external.get("appris", {}).get(transcript_base) if transcript_base else None,
                "gene_centric": accession == external.get("gene_centric_accession"),
                "duplicates": duplicates.get(accession, []),
                "features": analyse_features(record),
                "isoforms": extract_isoforms(record),
                "function_evidence": function_evidence(record),
                "links": external_links(record),
                "raw_record": record,
            }
        )
    return output


def _common_domain(records: list[dict[str, Any]]) -> dict[str, Any] | None:
    counts: Counter[str] = Counter()
    display: dict[str, str] = {}
    for record in records:
        unique = set()
        for domain in record.get("features", {}).get("domains", []) or []:
            name = str(domain.get("name") or "").strip()
            norm = normalize_text(name)
            if norm and norm not in unique:
                counts[norm] += 1
                display[norm] = name
                unique.add(norm)

    if not counts:
        return None
    ranked = counts.most_common()
    top_name, top_count = ranked[0]
    if len(ranked) > 1 and ranked[1][1] == top_count:
        return None
    if top_count < 2:
        return None
    return {"normalized": top_name, "name": display[top_name], "count": top_count}


def build_explanations(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not records:
        return []

    lengths = [record["length"] for record in records if record.get("length") is not None]
    min_length = min(lengths) if lengths else None
    max_length = max(lengths) if lengths else None
    existence_counts = Counter(record.get("existence") for record in records)
    name_counts = Counter(record.get("name") for record in records if record.get("name"))
    common_name = name_counts.most_common(1)[0][0] if name_counts else None
    signal_count = sum(bool(record.get("features", {}).get("signal")) for record in records)
    common_domain = _common_domain(records)

    explained: list[dict[str, Any]] = []
    for record in records:
        sentences: list[str] = []
        evidence: list[dict[str, str]] = []
        accession = record["accession"]

        if record.get("duplicates"):
            others = ", ".join(record["duplicates"])
            sentences.append(
                f"Its amino-acid sequence is exactly the same as {others}, even though UniProt lists them as separate accessions."
            )
            evidence.append({"claim": "Exact sequence match", "source": "Computed from UniProt sequences", "value": others})

        length = record.get("length")
        if length is not None and min_length is not None and max_length is not None and len(records) > 1:
            if min_length != max_length and length == min_length:
                sentences.append(f"At {length} amino acids, this is the shortest sequence among the {len(records)} records retrieved.")
            elif min_length != max_length and length == max_length:
                count = sum(item.get("length") == max_length for item in records)
                if count == 1:
                    sentences.append(f"At {length} amino acids, this is the longest sequence in the retrieved set.")
                else:
                    sentences.append(f"At {length} amino acids, this is one of {count} records with the longest sequence in the retrieved set.")

        if common_name and record.get("name") and record["name"] != common_name:
            sentences.append(f'UniProt names this entry "{record["name"]}", while the most common name in this retrieved set is "{common_name}".')

        if record.get("transcript"):
            if record.get("canonical"):
                sentence = f"UniProt links this protein to Ensembl transcript {record['transcript']}, and Ensembl marks that transcript as canonical."
                if record.get("appris"):
                    sentence += f" APPRIS labels it {record['appris']}."
                if record.get("gene_centric"):
                    sentence += " UniProt's gene-centric service also returns this accession as the representative entry for the gene group."
                sentences.append(sentence)
        else:
            sentences.append("UniProt does not return an Ensembl transcript cross-reference for this accession.")

        existence = record.get("existence")
        if len(records) > 1 and existence_counts.get(existence, 0) == 1:
            sentences.append(f"This is the only record in this set with UniProt Protein Existence classified as {existence}.")

        domains = record.get("features", {}).get("domains", []) or []
        domain_norms = {normalize_text(item.get("name")) for item in domains if item.get("name")}
        if common_domain:
            common_label = common_domain["name"]
            if common_domain["normalized"] in domain_norms:
                matching = next((item for item in domains if normalize_text(item.get("name")) == common_domain["normalized"]), None)
                pos = matching.get("position") if matching else None
                sentences.append(
                    f"UniProt annotates the {common_label} domain" + (f" at amino acids {pos}." if pos else ".")
                )
            else:
                sentences.append(
                    f"UniProt does not annotate the {common_label} domain on this sequence, while it is annotated on {common_domain['count']} of the {len(records)} retrieved records."
                )

        other_domains = [
            item for item in domains
            if not common_domain or normalize_text(item.get("name")) != common_domain["normalized"]
        ]
        if other_domains:
            descriptions = [
                str(item.get("name")) + (f" ({item.get('position')})" if item.get("position") else "")
                for item in other_domains
            ]
            sentences.append("Other UniProt domain annotations on this sequence: " + "; ".join(descriptions) + ".")

        signal = record.get("features", {}).get("signal", []) or []
        if signal:
            sentences.append(f"UniProt annotates a signal peptide at {signal[0]}.")
        elif 0 < signal_count < len(records):
            sentences.append(f"UniProt does not annotate a signal peptide on this sequence, while it does on {signal_count} of the {len(records)} retrieved records.")

        membrane = record.get("features", {}).get("transmembrane", []) or []
        if membrane:
            sentences.append(f"UniProt annotates a membrane-spanning region at amino acids {membrane[0]}.")

        disordered = [value for value in record.get("features", {}).get("disordered", []) or [] if value]
        if disordered:
            sentences.append(f"UniProt annotates {len(disordered)} disordered region{'s' if len(disordered) != 1 else ''}: {', '.join(disordered)}.")

        explained.append(
            {
                "accession": accession,
                "length": length,
                "sentences": sentences,
                "record": record,
                "evidence": evidence,
            }
        )

    return explained




def build_difference_summary(records: list[dict[str, Any]]) -> list[str]:
    """Summarise observable differences across a UniProt record set without ranking entries."""
    if not records or len(records) == 1:
        return []

    notes: list[str] = []
    lengths = sorted({record.get("length") for record in records if record.get("length") is not None})
    if len(lengths) > 1:
        notes.append(f"Sequence lengths vary from {lengths[0]} to {lengths[-1]} amino acids across the retrieved records.")

    names = {str(record.get("name") or "").strip() for record in records if record.get("name")}
    if len(names) > 1:
        notes.append(f"UniProt uses {len(names)} different protein-name labels across this record set.")

    transcript_count = sum(bool(record.get("transcript")) for record in records)
    if 0 < transcript_count < len(records):
        notes.append(f"{transcript_count} of {len(records)} records have an Ensembl transcript cross-reference; the others do not.")

    existence = {str(record.get("existence") or "").strip() for record in records if record.get("existence")}
    if len(existence) > 1:
        notes.append("Protein-existence evidence differs across the records (for example protein-level, transcript-level, or inferred evidence).")

    duplicate_pairs = sum(1 for record in records if record.get("duplicates"))
    if duplicate_pairs:
        notes.append(f"{duplicate_pairs} record{'s' if duplicate_pairs != 1 else ''} share an identical amino-acid sequence with at least one other accession in this set.")

    isoform_entries = sum(bool(record.get("isoforms")) for record in records)
    if isoform_entries:
        notes.append(f"{isoform_entries} entr{'ies' if isoform_entries != 1 else 'y'} describe one or more UniProt isoforms inside the entry itself.")

    return notes[:4]

def explain_records(
    records: list[dict[str, Any]],
    gene: str,
    species: str,
    client: DataClient,
    *,
    include_external: bool = False,
) -> dict[str, Any]:
    """Build the UniProt comparison without blocking on secondary databases.

    The Colab version queried Ensembl, APPRIS and UniProt gene-centric data as
    part of every explanation. On a serverless web request those secondary
    services can dominate latency. The core result now renders immediately
    from UniProt; secondary annotations can be requested separately.
    """
    warnings: list[dict[str, str]] = []
    external = {
        "ensembl_gene": None,
        "ensembl_canonical": None,
        "ensembl_species": None,
        "gene_centric_accession": None,
        "appris": {},
    }
    if include_external:
        external, warnings = get_external_annotation_context(records, client)

    comparison = build_comparison_records(records, external)
    explained = build_explanations(comparison)

    review_counts = Counter(item["review"]["status"] for item in comparison)
    reviewed = review_counts.get("reviewed", 0)
    unreviewed = review_counts.get("unreviewed", 0)
    total = len(comparison)

    review_summary = (
        f"{reviewed} of these {total} UniProtKB records are reviewed (Swiss-Prot); "
        f"{unreviewed} are unreviewed (TrEMBL)."
    )

    has_isoforms = any(item.get("isoforms") for item in comparison)
    difference_summary = build_difference_summary(comparison)

    return {
        "gene": gene,
        "species": species,
        "record_count": total,
        "reviewed_count": reviewed,
        "unreviewed_count": unreviewed,
        "review_summary": review_summary,
        "review_explanation": (
            "Reviewed (Swiss-Prot) entries are manually reviewed and annotated by UniProt curators. "
            "Unreviewed (TrEMBL) entries are computationally annotated and await full manual review."
        ),
        "isoform_explanation": ISOFORM_EXPLANATION if has_isoforms else None,
        "difference_summary": difference_summary,
        "records": explained,
        "warnings": warnings,
        "external_annotations": external,
        "external_annotations_loaded": include_external,
    }


def load_external_annotations(records: list[dict[str, Any]], client: DataClient) -> dict[str, Any]:
    external, warnings = get_external_annotation_context(records, client)
    return {"external_annotations": external, "warnings": warnings}
