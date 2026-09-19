# UniProt Explorer — v18 smooth web MVP

This build keeps the final Colab search/evidence rules as the scientific source of truth and improves the web product around them.

## What the user can do

- Search once in plain language: `TP53 in human`, `anthrax in humans`, `Li-Fraumeni syndrome`, etc.
- See direct protein/gene identity matches separately from evidence-backed discovery.
- Open the complete UniProt record set for a gene in its source organism.
- See reviewed (Swiss-Prot) vs unreviewed (TrEMBL), isoforms, evidence, domains and links.
- Filter a large record set by review status or text.
- Select records and download the selection.
- Download a whole gene record set or one entry as FASTA, TSV, JSON, XML or UniProt flat text.
- Open UniProt, AlphaFold DB and PDB links.
- Load Ensembl/APPRIS/gene-centric context only when wanted.

## Why v18 is faster than v17

The previous web adapter repeated work that the Colab notebook did not need to repeat in an interactive website. v18 keeps the search semantics but changes request timing:

1. A mentioned-organism UniProt search is performed once and reused for both direct-name matching and source-organism discovery.
2. If a direct identity match is found, the app immediately retrieves that gene's record set and does not wait for unrelated global/virus discovery.
3. If discovery is needed, virus-host and global branches run concurrently.
4. Ensembl, APPRIS and UniProt gene-centric calls are lazy and do not block the first result page.
5. Expensive spelling recovery is available on demand instead of delaying every zero-result search.
6. Repeated browser requests in the same session are cached client-side; HTTP responses are also cached server-side while the Vercel instance is warm.
7. Network retries are bounded so one slow upstream branch cannot hold the page for minutes.

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

Current suite: 38 tests.

## Deployment

The repository can be deployed directly to Vercel as a FastAPI project. `app/main.py` is the FastAPI entry point and also serves the static frontend.
