from __future__ import annotations

from typing import Any

from .explain import explain_records, load_external_annotations
from .http import DataClient, SourceUnavailable
from .query import concept_orthographic_variants, extract_search_concepts
from .records import get_gene_name, get_protein_name
from .search import (
    direct_search_phrases,
    discovery_search,
    fetch_all_records_for_gene,
    group_direct_matches,
    suggest_uniprot_wording_correction,
)
from .taxonomy import resolve_species, taxonomy_display_name
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


class ProteinService:
    """Web adapter around the final Colab search logic.

    The scientific/search behavior remains Colab-derived. The web adapter only
    changes *when* work is done so that a user is not forced to wait for
    duplicate searches or secondary databases before seeing UniProt results.
    """

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
        confirmed_spelling: str | None = None,
    ) -> dict[str, Any]:
        query = query.strip()
        if not query:
            return {"status": "error", "message": "Enter a protein, gene, or biological concept.", "warnings": []}

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
                    "warnings": _dedupe_warnings(warnings),
                }

            if species["status"] == "resolved":
                match = species["match"]
                record = match["record"]
                taxon_id = record.get("taxonId")
                detected_phrase = match.get("input_phrase")
                species_label = record.get("commonName") or record.get("scientificName")
                organism_resolution = "resolved"
            elif species["status"] == "broad_only":
                warnings.append(
                    {
                        "source": "organism resolver",
                        "message": "Part of the wording matched a broader taxonomic group, so the search continued without forcing one species.",
                    }
                )

        has_species = taxon_id is not None

        concept_info = extract_search_concepts(query, detected_phrase)
        original_strict_concept = normalize_text(str(concept_info.get("primary") or ""))
        strict_concept = normalize_text(confirmed_spelling or original_strict_concept)
        strict_variants = concept_orthographic_variants(strict_concept)
        related_levels = list(concept_info.get("related_levels") or [])

        if not strict_concept:
            return {
                "status": "error",
                "message": "I identified an organism but could not identify useful protein or biological concept wording.",
                "warnings": _dedupe_warnings(warnings),
            }

        def organism_payload() -> dict[str, Any] | None:
            if not has_species:
                return None
            return {
                "taxon_id": taxon_id,
                "name": species_label,
                "input_phrase": detected_phrase,
                "resolution": organism_resolution,
            }

        # ------------------------------------------------------------
        # STRICT SEARCH
        # ------------------------------------------------------------
        source_rows: list[dict[str, Any]] | None = None
        direct_groups: list[dict[str, Any]] = []

        if has_species:
            # One mentioned-organism search does double duty: it is used first
            # for direct-name matching and then reused for source discovery.
            direct = direct_search_phrases(strict_variants, self.client, taxon_id=taxon_id)
            warnings.extend(direct.get("warnings", []))
            source_rows = direct.get("search_records", [])
            direct_groups = group_direct_matches(direct.get("matches", []))

            # A direct identity result does not wait for virus/global discovery.
            if direct_groups:
                if len(direct_groups) == 1:
                    group = direct_groups[0]
                    return self.explain_gene(
                        gene=str(group["gene"]),
                        taxon_id=int(group["taxon_id"]),
                        species_name=str(group["organism"]),
                        inherited_warnings=warnings,
                        search_context={
                            "original_query": query,
                            "search_concept": strict_concept,
                            "mode": "strict",
                            "direct_reasons": group.get("reasons", []),
                        },
                    )

                return {
                    "status": "needs_protein_choice",
                    "message": (
                        f'I found more than one direct UniProt protein/gene-name match for "{strict_concept}" '
                        f'in {species_label or "the selected organism"}. Choose the protein you want to inspect.'
                    ),
                    "query": strict_concept,
                    "original_query": query,
                    "options": direct_groups,
                    "organism": organism_payload(),
                    "warnings": _dedupe_warnings(warnings),
                }

        strict_discovery = discovery_search(
            strict_variants,
            self.client,
            taxon_id=taxon_id if has_species else None,
            organism_name=species_label,
            source_rows=source_rows,
        )
        warnings.extend(strict_discovery.get("warnings", []))

        if strict_discovery.get("results"):
            return {
                "status": "no_direct_match",
                "query": strict_concept,
                "original_query": query,
                "message": (
                    f'I did not find a direct UniProt gene/protein-name match for "{strict_concept}"'
                    + (f" in {species_label}." if has_species and species_label else ".")
                    + " The results below are traceable discovery options, not a guess at one correct protein."
                ),
                "search_mode": "strict",
                "organism": organism_payload(),
                "discovery": strict_discovery["results"],
                "warnings": _dedupe_warnings(warnings),
            }

        # ------------------------------------------------------------
        # RELATED DISCOVERY
        # ------------------------------------------------------------
        for related_phrase in related_levels:
            related_variants = concept_orthographic_variants(related_phrase)
            related_source_rows = None
            if has_species:
                related_source = direct_search_phrases(related_variants, self.client, taxon_id=taxon_id)
                warnings.extend(related_source.get("warnings", []))
                related_source_rows = related_source.get("search_records", [])

            related = discovery_search(
                related_variants,
                self.client,
                taxon_id=taxon_id if has_species else None,
                organism_name=species_label,
                source_rows=related_source_rows,
            )
            warnings.extend(related.get("warnings", []))
            if related.get("results"):
                return {
                    "status": "no_direct_match",
                    "query": related_phrase,
                    "original_query": query,
                    "message": (
                        f'I did not find an exact or explainable UniProt match for the full wording '
                        f'"{original_strict_concept}". The results below use the broader related wording '
                        f'"{related_phrase}" and are not presented as exact matches for the original wording.'
                    ),
                    "search_mode": "related",
                    "related_from": original_strict_concept,
                    "organism": organism_payload(),
                    "discovery": related["results"],
                    "warnings": _dedupe_warnings(warnings),
                }

        # Spelling recovery still exists, but is no longer allowed to hold up
        # every zero-result web search. The UI can request it explicitly.
        return {
            "status": "no_direct_match",
            "query": strict_concept,
            "original_query": query,
            "message": (
                f'I did not find an exact or traceable UniProt match for "{original_strict_concept}". '
                "The tool did not invent an association from classification labels or unrelated metadata."
            ),
            "search_mode": "none",
            "organism": organism_payload(),
            "discovery": [],
            "can_check_spelling": bool(has_species and confirmed_spelling is None),
            "warnings": _dedupe_warnings(warnings),
        }

    def explain_gene(
        self,
        gene: str,
        taxon_id: int,
        *,
        species_name: str | None = None,
        inherited_warnings: list[dict[str, str]] | None = None,
        search_context: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        records, warnings = fetch_all_records_for_gene(gene, taxon_id, self.client)
        warnings = [*(inherited_warnings or []), *warnings]
        if not records:
            return {
                "status": "not_found",
                "gene": gene,
                "taxon_id": taxon_id,
                "message": "No UniProtKB records were returned for that gene and organism.",
                "warnings": _dedupe_warnings(warnings),
            }

        if species_name is None:
            organism = records[0].get("organism", {}) or {}
            species_name = organism.get("commonName") or organism.get("scientificName") or str(taxon_id)

        # Core explanations use UniProt only. Ensembl/APPRIS/gene-centric data
        # is available through a lazy endpoint and does not block first paint.
        payload = explain_records(records, gene, species_name, self.client, include_external=False)
        payload["status"] = "ready"
        payload["taxon_id"] = taxon_id
        payload["accessions"] = [
            str(record.get("primaryAccession"))
            for record in records
            if record.get("primaryAccession")
        ]
        payload["download_scope"] = {"mode": "gene", "gene": gene, "taxon_id": taxon_id}
        payload["warnings"] = _dedupe_warnings([*warnings, *payload.get("warnings", [])])
        if search_context:
            payload["search_context"] = search_context
        return payload

    def explain_accession(
        self,
        accession: str,
        *,
        species_name: str | None = None,
        inherited_warnings: list[dict[str, str]] | None = None,
    ) -> dict[str, Any]:
        warnings = list(inherited_warnings or [])
        try:
            record = self.client.uniprot_entry(accession)
        except SourceUnavailable as exc:
            warnings.append(exc.as_warning())
            return {
                "status": "not_found",
                "message": f"UniProt entry {accession} could not be retrieved.",
                "warnings": _dedupe_warnings(warnings),
            }

        organism = record.get("organism", {}) or {}
        species = species_name or organism.get("commonName") or organism.get("scientificName") or "the source organism"
        gene = get_gene_name(record) or get_protein_name(record) or accession
        payload = explain_records([record], str(gene), str(species), self.client, include_external=False)
        payload["status"] = "ready"
        payload["single_accession"] = accession
        payload["accessions"] = [accession]
        payload["download_scope"] = {"mode": "accessions", "accessions": [accession]}
        payload["warnings"] = _dedupe_warnings([*warnings, *payload.get("warnings", [])])
        return payload

    def discover(self, query: str, taxon_id: int | None = None) -> dict[str, Any]:
        query = normalize_text(query)
        if not query:
            return {"results": [], "warnings": [], "message": "Enter a biological concept."}
        payload = discovery_search(query, self.client, taxon_id=taxon_id)
        return {
            "results": payload["results"],
            "warnings": payload["warnings"],
            "query": query,
        }

    def check_spelling(self, phrase: str, taxon_id: int) -> dict[str, Any]:
        suggestion, warnings = suggest_uniprot_wording_correction(phrase, taxon_id, self.client)
        return {"suggestion": suggestion, "warnings": warnings}

    def external_context(self, gene: str, taxon_id: int) -> dict[str, Any]:
        records, warnings = fetch_all_records_for_gene(gene, taxon_id, self.client)
        if not records:
            return {
                "status": "not_found",
                "message": "No UniProtKB records were returned for that gene and organism.",
                "warnings": warnings,
            }
        payload = load_external_annotations(records, self.client)
        payload["status"] = "ready"
        payload["warnings"] = _dedupe_warnings([*warnings, *payload.get("warnings", [])])
        return payload
