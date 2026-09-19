from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Iterable

import requests
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from .core.http import UNIPROT_STREAM
from .core.service import ProteinService

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

app = FastAPI(title="UniProt Explorer", version="0.18.0")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=False,
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)

service = ProteinService()
STATIC = Path(__file__).parent / "static"
app.mount("/static", StaticFiles(directory=STATIC), name="static")


class LookupRequest(BaseModel):
    query: str = Field(min_length=1, max_length=300)
    taxon_id: int | None = None
    organism_name: str | None = None
    organism_phrase: str | None = None
    resolve_organism: bool = True
    confirmed_spelling: str | None = None


class ExplainRequest(BaseModel):
    gene: str = Field(min_length=1, max_length=100)
    taxon_id: int
    species_name: str | None = None


class EntryRequest(BaseModel):
    accession: str = Field(min_length=1, max_length=50)
    species_name: str | None = None


class DiscoveryRequest(BaseModel):
    query: str = Field(min_length=1, max_length=300)
    taxon_id: int | None = None


class SpellingRequest(BaseModel):
    phrase: str = Field(min_length=1, max_length=200)
    taxon_id: int


class ExternalRequest(BaseModel):
    gene: str = Field(min_length=1, max_length=100)
    taxon_id: int


class SelectedDownloadRequest(BaseModel):
    accessions: list[str] = Field(min_length=1, max_length=500)
    format: str = Field(default="fasta")


DOWNLOAD_FORMATS = {
    "fasta": ("fasta", "text/plain; charset=utf-8"),
    "tsv": ("tsv", "text/tab-separated-values; charset=utf-8"),
    "json": ("json", "application/json"),
    "xml": ("xml", "application/xml"),
    "txt": ("txt", "text/plain; charset=utf-8"),
}


def _safe_accession(value: str) -> str:
    value = str(value or "").strip()
    if not re.fullmatch(r"[A-Za-z0-9_.-]{1,50}", value):
        raise HTTPException(status_code=400, detail=f"Invalid UniProt accession: {value}")
    return value


def _download_response(query: str, fmt: str, filename: str):
    fmt = fmt.lower().strip()
    if fmt not in DOWNLOAD_FORMATS:
        raise HTTPException(status_code=400, detail="Format must be fasta, tsv, json, xml, or txt.")

    extension, media_type = DOWNLOAD_FORMATS[fmt]
    params = {"query": query, "format": fmt}
    if fmt == "tsv":
        params["fields"] = "accession,id,reviewed,protein_name,gene_names,organism_name,length,protein_existence"

    try:
        upstream = requests.get(
            UNIPROT_STREAM,
            params=params,
            timeout=45,
            stream=True,
            headers={"User-Agent": "UniProtExplorer/0.18", "Accept": "*/*"},
        )
        upstream.raise_for_status()
    except requests.RequestException as exc:
        logger.exception("UniProt download failed query=%s format=%s", query, fmt)
        raise HTTPException(status_code=502, detail=f"UniProt download was unavailable: {exc}") from exc

    def iterator() -> Iterable[bytes]:
        try:
            for chunk in upstream.iter_content(chunk_size=64 * 1024):
                if chunk:
                    yield chunk
        finally:
            upstream.close()

    headers = {"Content-Disposition": f'attachment; filename="{filename}.{extension}"'}
    return StreamingResponse(iterator(), media_type=media_type, headers=headers)


@app.exception_handler(Exception)
async def unexpected_error_handler(request, exc: Exception):
    logger.exception("Unhandled request failure for %s", request.url.path)
    return JSONResponse(
        status_code=500,
        content={
            "status": "error",
            "message": "The request could not be completed. Check the source notes or try again shortly.",
        },
    )


@app.get("/")
def home():
    return FileResponse(STATIC / "index.html")


@app.get("/api/health")
def health():
    return {"status": "ok", "version": "0.18.0"}


@app.post("/api/lookup")
def lookup(request: LookupRequest):
    return service.lookup(
        request.query,
        taxon_id=request.taxon_id,
        organism_name=request.organism_name,
        organism_phrase=request.organism_phrase,
        resolve_organism=request.resolve_organism,
        confirmed_spelling=request.confirmed_spelling,
    )


@app.post("/api/explain")
def explain(request: ExplainRequest):
    return service.explain_gene(
        request.gene,
        request.taxon_id,
        species_name=request.species_name,
    )


@app.post("/api/entry")
def entry(request: EntryRequest):
    return service.explain_accession(
        request.accession,
        species_name=request.species_name,
    )


@app.post("/api/discover")
def discover(request: DiscoveryRequest):
    return service.discover(request.query, taxon_id=request.taxon_id)


@app.post("/api/spelling")
def spelling(request: SpellingRequest):
    return service.check_spelling(request.phrase, request.taxon_id)


@app.post("/api/external-context")
def external_context(request: ExternalRequest):
    return service.external_context(request.gene, request.taxon_id)


@app.get("/api/download/gene")
def download_gene(
    gene: str = Query(min_length=1, max_length=100),
    taxon_id: int = Query(gt=0),
    format: str = Query(default="fasta"),
):
    safe_gene = re.sub(r"[^A-Za-z0-9_.-]+", "_", gene).strip("_") or "gene"
    query = f"(gene_exact:{gene}) AND (organism_id:{taxon_id})"
    return _download_response(query, format, f"uniprot_{safe_gene}_{taxon_id}")


@app.get("/api/download/entry/{accession}")
def download_entry(accession: str, format: str = Query(default="fasta")):
    accession = _safe_accession(accession)
    return _download_response(f"accession:{accession}", format, f"uniprot_{accession}")


@app.post("/api/download/selected")
def download_selected(request: SelectedDownloadRequest):
    accessions = list(dict.fromkeys(_safe_accession(value) for value in request.accessions))
    query = "(" + " OR ".join(f"accession:{value}" for value in accessions) + ")"
    return _download_response(query, request.format, f"uniprot_selected_{len(accessions)}")
