from __future__ import annotations

from typing import Any

from .explain import explain_records
from .http import DataClient
from .search import direct_search, discovery_search, fetch_all_records_for_gene, group_direct_matches
from .taxonomy import resolve_species, taxonomy_display_name
from .text import normalize_text


class ProteinService:
    def __init__(self, client: DataClient | None = None) -> None:
        self.client = client or DataClient()

    @staticmethod
    def _remove_species_phrase(text: str, phrase: str | None) -> str:
        if not phrase:
            return normalize_text(text)
        normalized = normalize_text(text)
        target = normalize_text(phrase)
        if not target:
            return normalized
        # Only remove a literal resolved phrase. No semantic rewriting.
        value = normalized.replace(target, " ")
        for word in (" in ", " from ", " of "):
            value = value.replace(word, " ")
        return " ".join(value.split())

    def lookup(
        self,
        query: str,
        *,
        taxon_id: int | None = None,
        organism_name: str | None = None,
        organism_phrase: str | None = None,
        resolve_organism: bool = True,
    ) -> dict[str, Any]:
        query = query.strip()
        if not query:
            return {"status": "error", "message": "Enter a protein or gene name.", "warnings": []}

        warnings: list[dict[str, str]] = []
        detected_phrase = None
        species_label = organism_name

        if taxon_id is None and resolve_organism:
            species = resolve_species(query, self.client)
            warnings.extend(species.get("warnings", []))
            if species["status"] == "ambiguous":
                return {
                    "status": "needs_organism_choice",
                    "message": "The wording could refer to more than one organism. Choose one, or continue without an organism.",
                    "options": [
                        {
                            "taxon_id": item["record"].get("taxonId"),
                            "label": taxonomy_display_name(item["record"]),
                            "rank": item["record"].get("rank"),
                            "input_phrase": item["input_phrase"],
                        }
                        for item in species["matches"]
                    ],
                    "warnings": warnings,
                }
            if species["status"] == "resolved":
                match = species["match"]
                record = match["record"]
                taxon_id = record.get("taxonId")
                detected_phrase = match.get("input_phrase")
                species_label = record.get("commonName") or record.get("scientificName")

        if organism_phrase:
            detected_phrase = organism_phrase

        target = self._remove_species_phrase(query, detected_phrase)
        if not target:
            return {
                "status": "error",
                "message": "I found an organism but not a protein or gene name to inspect.",
                "warnings": warnings,
            }

        direct = direct_search(target, self.client, taxon_id=taxon_id)
        warnings.extend(direct.get("warnings", []))
        groups = group_direct_matches(direct["matches"])

        if len(groups) == 1:
            group = groups[0]
            return self.explain_gene(
                gene=str(group["gene"]),
                taxon_id=int(group["taxon_id"]),
                species_name=str(group["organism"]),
                inherited_warnings=warnings,
            )

        if len(groups) > 1:
            return {
                "status": "needs_protein_choice",
                "message": "More than one direct UniProt gene/protein identity match was found. Choose the organism/gene you meant.",
                "options": groups,
                "warnings": warnings,
            }

        discovery = discovery_search(target, self.client, taxon_id=taxon_id)
        warnings.extend(discovery.get("warnings", []))
        return {
            "status": "no_direct_match",
            "query": target,
            "message": (
                "No direct UniProt gene/protein-name match was found. The related results below are discovery options, not a guess at one correct protein."
            ),
            "discovery": discovery["results"],
            "warnings": warnings,
        }

    def explain_gene(
        self,
        gene: str,
        taxon_id: int,
        *,
        species_name: str | None = None,
        inherited_warnings: list[dict[str, str]] | None = None,
    ) -> dict[str, Any]:
        records, warnings = fetch_all_records_for_gene(gene, taxon_id, self.client)
        warnings = [*(inherited_warnings or []), *warnings]
        if not records:
            return {
                "status": "not_found",
                "gene": gene,
                "taxon_id": taxon_id,
                "message": "No UniProtKB records were returned for that gene and organism.",
                "warnings": warnings,
            }

        if species_name is None:
            organism = records[0].get("organism", {}) or {}
            species_name = organism.get("commonName") or organism.get("scientificName") or str(taxon_id)

        payload = explain_records(records, gene, species_name, self.client)
        payload["status"] = "ready"
        payload["taxon_id"] = taxon_id
        payload["warnings"] = [*warnings, *payload.get("warnings", [])]
        return payload

    def discover(self, query: str, taxon_id: int | None = None) -> dict[str, Any]:
        payload = discovery_search(query, self.client, taxon_id=taxon_id)
        return {"status": "ready", **payload}
