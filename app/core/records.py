from __future__ import annotations

import json
import urllib.parse
from collections import defaultdict
from typing import Any

from .text import normalize_text


def get_gene_name(record: dict[str, Any]) -> str | None:
    genes = record.get("genes", []) or []
    if genes:
        gene_name = genes[0].get("geneName", {}) or {}
        if gene_name.get("value"):
            return str(gene_name["value"])
    return None


def get_all_gene_labels(record: dict[str, Any]) -> list[str]:
    output: list[str] = []
    for gene in record.get("genes", []) or []:
        for key in ("geneName", "orderedLocusName", "orfNames"):
            value = gene.get(key)
            if isinstance(value, dict) and value.get("value"):
                output.append(str(value["value"]))
            elif isinstance(value, list):
                for item in value:
                    if isinstance(item, dict) and item.get("value"):
                        output.append(str(item["value"]))
        for item in gene.get("synonyms", []) or []:
            if isinstance(item, dict) and item.get("value"):
                output.append(str(item["value"]))
    return list(dict.fromkeys(output))


def _name_value(node: Any) -> str | None:
    if isinstance(node, dict):
        if node.get("value"):
            return str(node["value"])
        for key in ("fullName", "shortName"):
            nested = node.get(key)
            if isinstance(nested, dict) and nested.get("value"):
                return str(nested["value"])
    return None


def get_protein_name(record: dict[str, Any]) -> str | None:
    description = record.get("proteinDescription", {}) or {}
    recommended = description.get("recommendedName", {}) or {}
    full = _name_value(recommended.get("fullName"))
    if full:
        return full
    submission = description.get("submissionNames", []) or []
    for item in submission:
        full = _name_value(item.get("fullName", {}) if isinstance(item, dict) else None)
        if full:
            return full
    return None


def get_alternative_names(record: dict[str, Any]) -> list[str]:
    output: list[str] = []
    description = record.get("proteinDescription", {}) or {}
    for item in description.get("alternativeNames", []) or []:
        if not isinstance(item, dict):
            continue
        for key in ("fullName", "shortNames"):
            node = item.get(key)
            if isinstance(node, dict):
                value = _name_value(node)
                if value:
                    output.append(value)
            elif isinstance(node, list):
                for child in node:
                    value = _name_value(child)
                    if value:
                        output.append(value)
    return list(dict.fromkeys(output))


def get_all_protein_search_names(record: dict[str, Any]) -> list[str]:
    names = [get_protein_name(record), *get_alternative_names(record)]
    return [name for name in dict.fromkeys(names) if name]


def direct_match_details(record: dict[str, Any], query: str) -> dict[str, Any] | None:
    """Only gene/protein identity fields can create a direct match."""
    target = normalize_text(query)
    if not target:
        return None

    gene_hits = [name for name in get_all_gene_labels(record) if normalize_text(name) == target]
    protein_hits = [name for name in get_all_protein_search_names(record) if normalize_text(name) == target]

    if not gene_hits and not protein_hits:
        return None
    return {
        "gene_matches": gene_hits,
        "protein_matches": protein_hits,
        "gene": get_gene_name(record),
        "protein_name": get_protein_name(record),
    }


def get_sequence(record: dict[str, Any]) -> str | None:
    sequence = record.get("sequence", {}) or {}
    value = sequence.get("value")
    return str(value) if value else None


def get_length(record: dict[str, Any]) -> int | None:
    sequence = record.get("sequence", {}) or {}
    value = sequence.get("length")
    if isinstance(value, int):
        return value
    seq = sequence.get("value")
    return len(seq) if seq else None


def get_protein_existence(record: dict[str, Any]) -> str | None:
    value = record.get("proteinExistence")
    return str(value) if value else None


def existence_kind(value: str | None) -> str:
    text = str(value or "").lower()
    if "protein level" in text:
        return "protein level"
    if "transcript level" in text:
        return "transcript level"
    if "homology" in text:
        return "inferred from homology"
    if "predicted" in text:
        return "predicted"
    if "uncertain" in text:
        return "uncertain"
    return str(value or "unknown")


def review_status(record: dict[str, Any]) -> dict[str, str]:
    entry_type = str(record.get("entryType") or "")
    lower = entry_type.lower()
    if "reviewed" in lower and "unreviewed" not in lower:
        return {
            "status": "reviewed",
            "label": "Reviewed (Swiss-Prot)",
            "explanation": "This UniProtKB entry has been manually reviewed and annotated by curators.",
        }
    if "unreviewed" in lower or "trembl" in lower:
        return {
            "status": "unreviewed",
            "label": "Unreviewed (TrEMBL)",
            "explanation": "This UniProtKB entry is computationally annotated and has not yet received full manual review.",
        }
    return {
        "status": "unknown",
        "label": entry_type or "Review status unavailable",
        "explanation": "UniProt did not expose a reviewed/unreviewed label for this record.",
    }


def get_proteome_ids(record: dict[str, Any]) -> list[str]:
    output: list[str] = []
    for xref in record.get("uniProtKBCrossReferences", []) or []:
        if xref.get("database") == "Proteomes" and xref.get("id"):
            output.append(str(xref["id"]))
    return output


def get_ensembl_details(record: dict[str, Any]) -> list[dict[str, str | None]]:
    rows: list[dict[str, str | None]] = []
    for xref in record.get("uniProtKBCrossReferences", []) or []:
        if xref.get("database") != "Ensembl":
            continue
        row: dict[str, str | None] = {"transcript": xref.get("id"), "protein": None, "gene": None}
        for prop in xref.get("properties", []) or []:
            key = prop.get("key")
            value = prop.get("value")
            if key == "ProteinId":
                row["protein"] = value
            elif key == "GeneId":
                row["gene"] = value
        rows.append(row)
    return rows


def feature_position(feature: dict[str, Any]) -> str | None:
    location = feature.get("location", {}) or {}
    start = (location.get("start", {}) or {}).get("value")
    end = (location.get("end", {}) or {}).get("value")
    if start is not None and end is not None:
        return f"{start}–{end}"
    return None


def analyse_features(record: dict[str, Any]) -> dict[str, list[Any]]:
    result: dict[str, list[Any]] = {
        "signal": [],
        "transmembrane": [],
        "domains": [],
        "disordered": [],
        "repeats": [],
    }
    for feature in record.get("features", []) or []:
        feature_type = feature.get("type")
        description = feature.get("description")
        position = feature_position(feature)
        if feature_type == "Signal":
            result["signal"].append(position)
        elif feature_type == "Transmembrane":
            result["transmembrane"].append(position)
        elif feature_type == "Domain":
            result["domains"].append({"name": description or "Unnamed domain", "position": position})
        elif feature_type == "Repeat":
            result["repeats"].append({"name": description or "Repeat", "position": position})
        elif feature_type == "Region" and description and "disorder" in description.lower():
            result["disordered"].append(position)
    return result


def extract_isoforms(record: dict[str, Any]) -> list[dict[str, Any]]:
    for comment in record.get("comments", []) or []:
        if comment.get("commentType") != "ALTERNATIVE PRODUCTS":
            continue
        output: list[dict[str, Any]] = []
        for isoform in comment.get("isoforms", []) or []:
            name = (isoform.get("name", {}) or {}).get("value")
            synonyms = [
                item.get("value")
                for item in isoform.get("synonyms", []) or []
                if isinstance(item, dict) and item.get("value")
            ]
            output.append(
                {
                    "name": name,
                    "ids": isoform.get("isoformIds", []) or [],
                    "synonyms": synonyms,
                    "sequence_status": isoform.get("isoformSequenceStatus"),
                }
            )
        return output
    return []


def pdb_ids(record: dict[str, Any]) -> list[str]:
    ids = [
        str(xref.get("id"))
        for xref in record.get("uniProtKBCrossReferences", []) or []
        if xref.get("database") == "PDB" and xref.get("id")
    ]
    return list(dict.fromkeys(ids))


def external_links(record: dict[str, Any]) -> dict[str, Any]:
    accession = str(record.get("primaryAccession") or "")
    pdb = pdb_ids(record)
    rcsb_query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": accession},
        },
        "return_type": "entry",
    }
    return {
        "uniprot": f"https://www.uniprot.org/uniprotkb/{accession}/entry" if accession else None,
        "alphafold": f"https://alphafold.ebi.ac.uk/entry/{accession}" if accession else None,
        "pdb": [
            {"id": pdb_id, "url": f"https://www.rcsb.org/structure/{pdb_id}"}
            for pdb_id in pdb
        ],
        "pdb_search": (
            "https://www.rcsb.org/search?request=" + urllib.parse.quote(json.dumps(rcsb_query, separators=(",", ":")))
            if accession
            else None
        ),
    }


def group_sequence_duplicates(records: list[dict[str, Any]]) -> dict[str, list[str]]:
    groups: dict[str, list[str]] = defaultdict(list)
    for record in records:
        seq = get_sequence(record)
        accession = record.get("primaryAccession")
        if seq and accession:
            groups[seq].append(str(accession))

    output: dict[str, list[str]] = {}
    for accessions in groups.values():
        if len(accessions) > 1:
            for accession in accessions:
                output[accession] = [item for item in accessions if item != accession]
    return output
