from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import timedelta
from pathlib import Path
from typing import Any, Iterable

import requests
import requests_cache
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger(__name__)

UNIPROT_TAXONOMY = "https://rest.uniprot.org/taxonomy/search"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb"
UNIPROT_GENECENTRIC = "https://www.ebi.ac.uk/proteins/api/genecentric"
ENSEMBL_LOOKUP = "https://rest.ensembl.org/lookup/id"
APPRIS_EXPORTER = "https://apprisws.bioinfo.cnio.es/rest/exporter/id"


class SourceUnavailable(RuntimeError):
    def __init__(self, source: str, message: str, *, url: str | None = None):
        super().__init__(message)
        self.source = source
        self.url = url
        self.message = message

    def as_warning(self) -> dict[str, str]:
        payload = {
            "source": self.source,
            "message": f"{self.source} was unavailable: {self.message}",
        }
        if self.url:
            payload["url"] = self.url
        return payload


class DataClient:
    """HTTP access layer with persistent caching, retries, and explicit failures."""

    def __init__(
        self,
        cache_name: str | None = None,
        cache_days: int = 7,
        timeout: int = 30,
        max_workers: int = 6,
    ) -> None:
        if cache_name is None:
            cache_dir = Path(".cache")
            cache_dir.mkdir(exist_ok=True)
            cache_name = str(cache_dir / "uniprot_http")

        self.timeout = timeout
        self.max_workers = max_workers
        self.session = requests_cache.CachedSession(
            cache_name=cache_name,
            backend="sqlite",
            expire_after=timedelta(days=cache_days),
            allowable_methods=("GET",),
            stale_if_error=True,
        )

        retry = Retry(
            total=3,
            connect=3,
            read=3,
            status=3,
            backoff_factor=0.6,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset({"GET"}),
            respect_retry_after_header=True,
        )
        adapter = HTTPAdapter(max_retries=retry, pool_connections=20, pool_maxsize=20)
        self.session.mount("https://", adapter)
        self.session.mount("http://", adapter)
        self.session.headers.update(
            {
                "User-Agent": "UniProtBeginnerExplainer/0.15 (educational prototype)",
                "Accept": "application/json",
            }
        )

    def _get_json(
        self,
        source: str,
        url: str,
        *,
        params: dict[str, Any] | None = None,
        headers: dict[str, str] | None = None,
    ) -> Any:
        try:
            response = self.session.get(
                url,
                params=params,
                headers=headers,
                timeout=self.timeout,
            )
            response.raise_for_status()
            return response.json()
        except (requests.RequestException, ValueError) as exc:
            logger.exception("%s request failed: %s params=%s", source, url, params)
            raise SourceUnavailable(source, str(exc), url=url) from exc

    def taxonomy_search(self, term: str, size: int = 25) -> list[dict[str, Any]]:
        data = self._get_json(
            "UniProt Taxonomy",
            UNIPROT_TAXONOMY,
            params={"query": term, "format": "json", "size": size},
        )
        return data.get("results", []) if isinstance(data, dict) else []

    def uniprot_search(self, query: str, size: int = 100) -> list[dict[str, Any]]:
        data = self._get_json(
            "UniProtKB search",
            UNIPROT_SEARCH,
            params={"query": query, "format": "json", "size": size},
        )
        return data.get("results", []) if isinstance(data, dict) else []

    def uniprot_stream(self, query: str) -> list[dict[str, Any]]:
        """Retrieve full UniProtKB records in one streamed query instead of N entry calls."""
        data = self._get_json(
            "UniProtKB stream",
            UNIPROT_STREAM,
            params={"query": query, "format": "json"},
        )
        return data.get("results", []) if isinstance(data, dict) else []

    def uniprot_entry(self, accession: str) -> dict[str, Any]:
        return self._get_json(
            "UniProtKB entry",
            f"{UNIPROT_ENTRY}/{accession}.json",
        )

    def gene_centric(self, proteome_id: str, gene_id: str) -> Any:
        return self._get_json(
            "UniProt gene-centric",
            UNIPROT_GENECENTRIC,
            params={"proteome": proteome_id, "gene": gene_id},
        )

    def ensembl_lookup(self, ensembl_gene: str) -> dict[str, Any]:
        return self._get_json(
            "Ensembl",
            f"{ENSEMBL_LOOKUP}/{ensembl_gene}",
            headers={"Accept": "application/json"},
        )

    def appris(self, species: str, ensembl_gene: str) -> Any:
        return self._get_json(
            "APPRIS",
            f"{APPRIS_EXPORTER}/{species}/{ensembl_gene}",
            params={"format": "json", "sc": "ensembl", "methods": "appris"},
        )

    def concurrent_taxonomy_search(
        self, terms: Iterable[str], size: int = 25
    ) -> tuple[dict[str, list[dict[str, Any]]], list[dict[str, str]]]:
        return self._concurrent(
            terms,
            lambda term: self.taxonomy_search(term, size=size),
        )

    def concurrent_uniprot_search(
        self, queries: Iterable[str], size: int = 100
    ) -> tuple[dict[str, list[dict[str, Any]]], list[dict[str, str]]]:
        return self._concurrent(
            queries,
            lambda query: self.uniprot_search(query, size=size),
        )

    def _concurrent(self, values: Iterable[str], fn):
        unique = list(dict.fromkeys(value for value in values if value))
        results: dict[str, Any] = {}
        warnings: list[dict[str, str]] = []

        with ThreadPoolExecutor(max_workers=min(self.max_workers, max(1, len(unique)))) as pool:
            futures = {pool.submit(fn, value): value for value in unique}
            for future in as_completed(futures):
                value = futures[future]
                try:
                    results[value] = future.result()
                except SourceUnavailable as exc:
                    warnings.append(exc.as_warning())
                    results[value] = []
                except Exception as exc:  # defensive: never silently swallow unknown failures
                    logger.exception("Unexpected concurrent source failure for %s", value)
                    warnings.append(
                        {
                            "source": "internal search worker",
                            "message": f"A search branch failed for {value!r}: {exc}",
                        }
                    )
                    results[value] = []

        return results, warnings
