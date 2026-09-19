# UniProt Beginner Explainer — v15 web refactor

This is the web-ready refactor of the single-cell UniProt prototype.

## What changed

- **No `input()` or `print()` in the core.** Core functions take arguments and return dictionaries.
- **Three layers:** `app/core` (library), `app/main.py` (FastAPI), `app/static` (thin frontend).
- **Persistent HTTP cache:** `requests-cache`, seven-day TTL by default.
- **Concurrent independent searches:** taxonomy and UniProt search variants run concurrently.
- **Batch full-record retrieval:** gene record sets use UniProtKB's `/stream` endpoint instead of fetching each accession serially.
- **No runtime NLTK install/download.** Organism resolution uses UniProt Taxonomy scientific/common/synonym names only. If it cannot resolve confidently, it does not guess.
- **Visible source failures:** source errors are logged and returned in `warnings` rather than silently shortening results.
- **Reviewed vs unreviewed:** every result exposes Swiss-Prot/TrEMBL status and the result set starts with a plain-English count.
- **Isoforms:** alternative products are exposed and explained separately from separate accessions.
- **Evidence codes:** UniProt FUNCTION evidence is surfaced, with conservative plain-English categories for experimental, similarity-based, computational/model, and imported evidence.
- **Generic domain comparison:** the old hard-coded VWFA comparison is gone. The explainer finds a genuinely common domain in the retrieved set and compares records against it only when the comparison is meaningful.
- **Prominent links:** UniProt, AlphaFold DB, and PDB/PDB search links appear on each record.
- **Golden tests:** 20+ tests enforce the provenance/no-invention behavior and the new beginner-facing features.

## Run locally

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
uvicorn app.main:app --reload
```

Then open `http://127.0.0.1:8000`.

## Run tests

```bash
pytest -q
```

## API

- `POST /api/lookup` — front-door protein/gene lookup.
- `POST /api/explain` — explain all UniProt records for a known gene + taxon.
- `POST /api/discover` — separate evidence-traceable discovery mode.
- `GET /api/health` — health check.

The API deliberately returns choice states instead of prompting in the middle of execution. For example, ambiguous organism wording returns `status: needs_organism_choice` with options for the frontend to render.

## Important behavior retained from the prototype

A direct match is still restricted to UniProt gene/protein-name fields. Function, disease, pathway, literature, or linked-database text cannot silently become a direct identity match. Discovery results are separate and include the exact UniProt-exposed field that caused them to appear.
