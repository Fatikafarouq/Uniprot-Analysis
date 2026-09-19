# UniProt Explorer — v20 consolidated MVP

This is the consolidated web build. It keeps the final Colab search/evidence rules as the scientific source of truth, while the web layer focuses on speed, clarity, downloads, and practical record comparison.

## Core behavior

- One plain-language search box: `TP53 in human`, `anthrax in humans`, `Li-Fraumeni syndrome`, etc.
- Conservative organism resolution using UniProt Taxonomy.
- Exact direct identity is preferred over incidental longer-name mentions.
- If there is no exact identity, the Colab contained-name fallback and evidence-backed discovery remain available.
- Direct identity, evidence-backed discovery, and unsupported/no-evidence results remain clearly separated.
- A discovery result with a gene can open the complete UniProt record set for that gene in its source organism.
- Evidence from separate fields is never combined to invent a biological relationship.

## Record explanation

The first record screen now foregrounds the explanation rather than database clutter:

- total records;
- reviewed vs unreviewed counts;
- why the query was treated as a direct identity match;
- a concise summary of what differs across the record set;
- per-record comparison notes;
- isoform information;
- function evidence/ECO interpretation;
- UniProt, AlphaFold and compact PDB access.

Large PDB lists are deliberately collapsed to one `PDB structures (N)` link instead of rendering hundreds of structure links in every card.

## Downloads

Users can download:

- the complete gene record set;
- one individual UniProt entry;
- a user-selected subset of records;

in FASTA, TSV, JSON, XML, or UniProt flat text.

## Performance safeguards

- direct identity searches stop before unrelated discovery work;
- mentioned-organism search rows are reused for discovery instead of fetched twice;
- independent discovery branches run concurrently;
- Ensembl/APPRIS/gene-centric context is lazy;
- HTTP responses are cached;
- retries are bounded to one retry for transient upstream failures;
- normal API calls use a 10-second upstream timeout;
- the browser stops a search after 45 seconds instead of leaving the user waiting for minutes;
- repeated requests in the same browser session are cached in memory.

## Project layout

```text
app/
  core/
    evidence.py
    explain.py
    http.py
    query.py
    records.py
    search.py
    service.py
    taxonomy.py
    text.py
  static/
    index.html
    app.js
    styles.css
  main.py
tests/
  test_golden.py
requirements.txt
```

## Tests

Run:

```bash
PYTHONPATH=. pytest -q
```

Current suite: **43 tests**.

The suite includes regression coverage for false evidence combination, literature provenance, reviewed/unreviewed classification, isoforms, ECO interpretation, generic domain comparison, exact-vs-contained identity matching, downloads, retry/cache recovery, record-set difference summaries, and compact UI handling for large PDB sets.

## Deployment

The repository can be deployed directly to Vercel as a FastAPI project. `app/main.py` is the FastAPI entry point and also serves the static frontend.
