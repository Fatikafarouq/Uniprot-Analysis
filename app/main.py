from __future__ import annotations

import logging
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from .core.service import ProteinService

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

app = FastAPI(title="UniProt Beginner Explainer", version="0.16.0")
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


class ExplainRequest(BaseModel):
    gene: str = Field(min_length=1, max_length=100)
    taxon_id: int
    species_name: str | None = None


class DiscoveryRequest(BaseModel):
    query: str = Field(min_length=1, max_length=300)
    taxon_id: int | None = None


@app.get("/")
def home():
    return FileResponse(STATIC / "index.html")


@app.get("/api/health")
def health():
    return {"status": "ok"}


@app.post("/api/lookup")
def lookup(request: LookupRequest):
    return service.lookup(
        request.query,
        taxon_id=request.taxon_id,
        organism_name=request.organism_name,
        organism_phrase=request.organism_phrase,
        resolve_organism=request.resolve_organism,
    )


@app.post("/api/explain")
def explain(request: ExplainRequest):
    return service.explain_gene(
        request.gene,
        request.taxon_id,
        species_name=request.species_name,
    )


@app.post("/api/discover")
def discover(request: DiscoveryRequest):
    return service.discover(request.query, taxon_id=request.taxon_id)
