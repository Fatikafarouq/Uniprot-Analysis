# UniProt Beginner Explainer — v17 Colab-parity web build

This build treats the final Colab prototype (`uniprot_single_cell_v16_explained_discovery`) as the source of truth for search/discovery behavior, then layers the web/product improvements on top.

## Core behavior preserved from Colab

- A single query can contain a protein/gene, organism, disease, phenotype, pathway, or conversational wording.
- Organism resolution is conservative. A taxonomy hit validates a candidate; it does not silently invent organism intent.
- The confirmed organism phrase is removed while preserving the biological concept.
- Strict concept wording is searched first. Broader right-hand wording is used only after the strict concept fails and is explicitly labelled as related discovery.
- Direct identity matches come only from UniProt protein/gene name fields. Disease/function/pathway/literature evidence cannot silently become identity evidence.
- In a resolved organism, direct protein-name matches are grouped by gene. Selecting a gene fetches **all UniProtKB records for that gene + organism**.
- Discovery keeps the three Colab paths: proteins from the mentioned organism, viral proteins linked by the host field, and global UniProt connections.
- A discovery result remains selectable. If it has a gene + source taxon, selection opens **all UniProt records for that gene in that source organism**; otherwise it opens the individual accession.
- Every discovery result explains exactly why it appeared; the tool does not invent a biological connection from unrelated metadata.

## Product/engineering improvements layered on top

- Core library + FastAPI API + frontend; no `input()` or `print()` in request logic.
- `requests-cache` with a seven-day TTL locally; Vercel uses an in-memory cache because the deployed filesystem is read-only. The cache survives warm instances but is not a guaranteed cross-instance persistent cache.
- Concurrent independent taxonomy/search requests.
- Full gene record sets use UniProtKB `/stream` rather than serial accession fetches.
- No runtime NLTK install/download. Organism resolution uses UniProt Taxonomy names/synonyms; spelling recovery is UniProt-grounded.
- Source failures are logged and exposed to the user as warnings.
- Reviewed (Swiss-Prot) vs unreviewed (TrEMBL) counts and labels.
- Plain-English isoform explanation and exposed UniProt alternative products.
- ECO evidence on FUNCTION annotations is surfaced conservatively.
- Generic common-domain comparison; no VWFA-specific branch.
- UniProt, AlphaFold DB, and PDB links on record cards.
- Golden tests protect no-invention rules and Colab/web parity.

## API

- `POST /api/lookup` — single front door; organism resolution, direct identity, or evidence-backed discovery.
- `POST /api/explain` — fetch and explain all UniProt records for a gene + taxon.
- `POST /api/entry` — inspect one accession when no gene-level expansion is possible.
- `POST /api/discover` — explicit discovery-only mode.
- `GET /api/health` — health/version check.

## Tests

```bash
pytest -q
```

The current suite contains 32 golden/parity tests.
