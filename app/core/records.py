from __future__ import annotations

import json
import urllib.parse
from collections import defaultdict
from typing import Any

from .text import normalize_text


DIRECT_TOKEN_STOPWORDS = {
    "a", "an", "and", "the", "of", "in", "for", "from",
    "with", "to", "on", "by", "protein", "proteins", "gene", "genes",
}


def get_gene_name(record: dict[str, Any]) -> str | None:
    genes = record.get("genes", []) or []
    if genes:
        gene_name = genes[0].get("geneName", {}) or {}
        if gene_name.get("value"):
            return str(gene_name["value"])
    return None


def get_all_gene_labels(record: dict[str, Any]) -> list[str]:
    """Mirror the Colab v14 identity fields exactly: names, synonyms, locus and ORF names."""
    labels: list[str] = []
    for gene_block in record.get("genes", []) or []:
        for key in ("geneName", "synonyms", "orderedLocusNames", "orfNames"):
            value = gene_block.get(key)
            if isinstance(value, dict):
                value = [value]
            if not isinstance(value, list):
                continue
            for item in value:
                label = item.get("value") if isinstance(item, dict) else item
                if label:
                    labels.append(str(label))
    return list(dict.fromkeys(labels))


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
    for item in description.get("submissionNames", []) or []:
        if not isinstance(item, dict):
            continue
        full = _name_value(item.get("fullName"))
        if full:
            return full
    return None


def _add_name_block(names: list[str], block: Any) -> None:
    if not isinstance(block, dict):
        return
    full = block.get("fullName")
    if isinstance(full, dict) and full.get("value"):
        names.append(str(full["value"]))
    for short in block.get("shortNames", []) or []:
        if isinstance(short, dict) and short.get("value"):
            names.append(str(short["value"]))


def get_all_protein_search_names(record: dict[str, Any]) -> list[str]:
    """Port of the Colab v14 protein-name collector, including components and short names."""
    description = record.get("proteinDescription", {}) or {}
    names: list[str] = []

    _add_name_block(names, description.get("recommendedName"))
    for key in ("submissionNames", "alternativeNames"):
        for block in description.get(key, []) or []:
            _add_name_block(names, block)

    for component_key in ("includes", "contains"):
        for component in description.get(component_key, []) or []:
            if not isinstance(component, dict):
                continue
            _add_name_block(names, component.get("recommendedName"))
            for key in ("submissionNames", "alternativeNames"):
                for block in component.get(key, []) or []:
                    _add_name_block(names, block)

    return list(dict.fromkeys(names))


def get_alternative_names(record: dict[str, Any]) -> list[str]:
    output: list[str] = []
    description = record.get("proteinDescription", {}) or {}
    for item in description.get("alternativeNames", []) or []:
        if not isinstance(item, dict):
            continue
        full = item.get("fullName")
        if isinstance(full, dict) and full.get("value"):
            output.append(str(full["value"]))
        for short in item.get("shortNames", []) or []:
            if isinstance(short, dict) and short.get("value"):
                output.append(str(short["value"]))
    return list(dict.fromkeys(output))


def meaningful_query_tokens(text: str) -> list[str]:
    return [
        token
        for token in normalize_text(text).split()
        if token not in DIRECT_TOKEN_STOPWORDS and len(token) > 1
    ]


def direct_match_details(record: dict[str, Any], query: str) -> dict[str, Any] | None:
    """Colab-v14 identity rule: only name fields qualify, with conservative name-only relaxation."""
    target = normalize_text(query)
    if not target:
        return None

    gene_labels = get_all_gene_labels(record)
    protein_names = get_all_protein_search_names(record)
    reasons: list[str] = []
    gene_hits: list[str] = []
    protein_hits: list[str] = []

    for label in gene_labels:
        if normalize_text(label) == target:
            gene_hits.append(label)
            reasons.append(f'Gene name matches "{label}" exactly.')

    for name in protein_names:
        if normalize_text(name) == target:
            protein_hits.append(name)
            reasons.append(f'Protein name matches "{name}" exactly.')

    if not reasons and len(target) >= 3:
        for name in protein_names:
            candidate = normalize_text(name)
            if target in candidate:
                protein_hits.append(name)
                reasons.append(
                    f'The wording appears directly in the UniProt protein name "{name}".'
                )

    if not reasons:
        query_tokens = meaningful_query_tokens(target)
        if len(set(query_tokens)) >= 2:
            query_set = set(query_tokens)
            for name in protein_names:
                candidate_tokens = set(meaningful_query_tokens(name))
                if query_set.issubset(candidate_tokens):
                    protein_hits.append(name)
                    reasons.append(
                        f'The key words occur in the UniProt protein name "{name}".'
                    )

    if not reasons:
        return None

    return {
        "gene_matches": list(dict.fromkeys(gene_hits)),
        "protein_matches": list(dict.fromkeys(protein_hits)),
        "gene": get_gene_name(record),
        "protein_name": get_protein_name(record),
        "reasons": list(dict.fromkeys(reasons)),
    }


def get_function_notes(record: dict[str, Any]) -> list[str]:
    output: list[str] = []
    for comment in record.get("comments", []) or []:
        if comment.get("commentType") != "FUNCTION":
            continue
        for text in comment.get("texts", []) or []:
            if isinstance(text, dict) and text.get("value"):
                output.append(str(text["value"]))
    return list(dict.fromkeys(output))


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
