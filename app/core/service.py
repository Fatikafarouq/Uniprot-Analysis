from __future__ import annotations

from typing import Any

from .explain import explain_records
from .http import DataClient
from .query import concept_orthographic_variants, extract_search_concepts
from .search import (
    direct_search_phrases,
    discovery_search,
    fetch_all_records_for_gene,
    group_direct_matches,
)
from .taxonomy import resolve_species, taxonomy_display_name
from .text import normalize_text


class ProteinService:
    def __init__(self, client: DataClient | None = None) -> None:
        self.client = client or DataClient()

    def lookup(
        self,
        query: str,
        *,
        taxon_id: int | None = None,
        organism_name: str | None = None,
        organism_phrase: str | None = None,
        resolve_organism: bool = True,
    ) -> dict[str, Any]:
        """One front door: resolve organism -> direct identity -> evidence-backed discovery."""
        query = query.strip()
        if not query:
            return {
                "status": "error",
                "message": "Enter a protein, gene, or biological concept.",
                "warnings": [],
            }

        warnings: list[dict[str, str]] = []
        detected_phrase = organism_phrase
        species_label = organism_name
        organism_resolution = "provided" if taxon_id is not None else None

        if taxon_id is None and resolve_organism:
            species = resolve_species(query, self.client)
            warnings.extend(species.get("warnings", []))

            if species["status"] == "ambiguous":
                return {
                    "status": "needs_organism_choice",
                    "message": (
                        "The wording could refer to more than one organism. "
                        "Choose one, or continue without an organism."
                    ),
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
                organism_resolution = "resolved"
            elif species["status"] == "broad_only":
                # A broad taxonomic hit is not silently turned into a species.
                warnings.append(
                    {
                        "source": "organism resolver",
                        "message": (
                            "Part of the wording matched a broader taxonomic group, "
                            "so the search continued without forcing one species."
                        ),
                    }
                )

        concept_info = extract_search_concepts(query, detected_phrase)
        strict_concept = str(concept_info["primary"] or "")
        strict_variants = list(concept_info["strict_variants"] or [])
        related_levels = list(concept_info["related_levels"] or [])

        if not strict_concept:
            return {
                "status": "error",
                "message": (
                    "I identified an organism but could not identify useful "
                    "protein or biological concept wording."
                ),
                "warnings": warnings,
            }

        # ------------------------------------------------------------
        # 1. DIRECT IDENTITY: only exact gene/protein name fields qualify.
        # ------------------------------------------------------------
        direct = direct_search_phrases(strict_variants, self.client, taxon_id=taxon_id)
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
                "message": (
                    "More than one direct UniProt gene/protein identity match was found. "
                    "Choose the organism/gene you meant."
                ),
                "options": groups,
                "warnings": warnings,
            }

        # ------------------------------------------------------------
        # 2. STRICT DISCOVERY: same biological wording, multiple transparent
        #    paths. If an organism was mentioned, do NOT restrict discovery to
        #    that organism only: also inspect virus-host and global UniProt.
        # ------------------------------------------------------------
        discovery = discovery_search(
            strict_variants,
            self.client,
            taxon_id=taxon_id,
            organism_name=species_label,
        )
        warnings.extend(discovery.get("warnings", []))

        if discovery["results"]:
            return {
                "status": "no_direct_match",
                "query": strict_concept,
                "original_query": query,
                "message": (
                    "No direct UniProt gene/protein-name match was found. "
                    "The results below are traceable discovery options, not a guess at one correct protein."
                ),
                "search_mode": "strict",
                "organism": {
                    "taxon_id": taxon_id,
                    "name": species_label,
                    "input_phrase": detected_phrase,
                    "resolution": organism_resolution,
                }
                if taxon_id is not None
                else None,
                "discovery": discovery["results"],
                "warnings": _dedupe_warnings(warnings),
            }

        # ------------------------------------------------------------
        # 3. RELATED DISCOVERY: only after the full concept fails. Shorter
        #    right-hand phrases are explicitly labelled as broader, never as
        #    exact replacements for what the user asked.
        # ------------------------------------------------------------
        for related_phrase in related_levels:
            related_variants = concept_orthographic_variants(related_phrase)
            related = discovery_search(
                related_variants,
                self.client,
                taxon_id=taxon_id,
                organism_name=species_label,
            )
            warnings.extend(related.get("warnings", []))
            if related["results"]:
                return {
                    "status": "no_direct_match",
                    "query": related_phrase,
                    "original_query": query,
                    "message": (
                        f'I did not find an exact or explainable UniProt match for the full wording '
                        f'"{strict_concept}". The results below use the broader related wording '
                        f'"{related_phrase}" and are labelled as discovery options.'
                    ),
                    "search_mode": "related",
                    "related_from": strict_concept,
                    "organism": {
                        "taxon_id": taxon_id,
                        "name": species_label,
                        "input_phrase": detected_phrase,
                        "resolution": organism_resolution,
                    }
                    if taxon_id is not None
                    else None,
                    "discovery": related["results"],
                    "warnings": _dedupe_warnings(warnings),
                }

        return {
            "status": "no_direct_match",
            "query": strict_concept,
            "original_query": query,
            "message": (
                f'I did not find an exact or traceable UniProt match for "{strict_concept}". '
                "The tool did not invent an association from unrelated metadata."
            ),
            "search_mode": "none",
            "organism": {
                "taxon_id": taxon_id,
                "name": species_label,
                "input_phrase": detected_phrase,
                "resolution": organism_resolution,
            }
            if taxon_id is not None
            else None,
            "discovery": [],
            "warnings": _dedupe_warnings(warnings),
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
        # Secondary discovery-only screen remains available, but the main box
        # now reaches the same evidence-backed discovery automatically.
        concept = normalize_text(query)
        payload = discovery_search(concept, self.client, taxon_id=taxon_id)
        return {"status": "ready", **payload}


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
