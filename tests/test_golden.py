from __future__ import annotations

from app.core.evidence import annotation_match_contexts, evidence_label, function_evidence
from app.core.explain import ISOFORM_EXPLANATION, build_comparison_records, build_explanations
from app.core.records import (
    direct_match_details,
    external_links,
    extract_isoforms,
    group_sequence_duplicates,
    review_status,
)
from app.core.search import fetch_all_records_for_gene
from app.core.taxonomy import candidate_species_phrases, taxonomy_record_is_placeholder


def base_record(accession="P1", gene="ABC1", name="Alpha protein", entry_type="UniProtKB reviewed (Swiss-Prot)"):
    return {
        "entryType": entry_type,
        "primaryAccession": accession,
        "proteinExistence": "1: Evidence at protein level",
        "organism": {"scientificName": "Homo sapiens", "commonName": "Human", "taxonId": 9606},
        "proteinDescription": {"recommendedName": {"fullName": {"value": name}}},
        "genes": [{"geneName": {"value": gene}, "synonyms": [{"value": f"{gene}S"}]}],
        "sequence": {"length": 100, "value": "A" * 100},
        "comments": [],
        "features": [],
        "uniProtKBCrossReferences": [],
        "references": [],
        "keywords": [],
    }


class StreamClient:
    def __init__(self, exact=None, fallback=None):
        self.exact = exact or []
        self.fallback = fallback or []
        self.calls = []

    def uniprot_stream(self, query):
        self.calls.append(query)
        return self.exact if "gene_exact:" in query else self.fallback


# 1
def test_no_direct_match_from_function_text():
    record = base_record()
    record["comments"] = [{"commentType": "FUNCTION", "texts": [{"value": "Important in DNA repair"}]}]
    assert direct_match_details(record, "DNA repair") is None


# 2
def test_direct_gene_name_match():
    record = base_record(gene="TP53")
    match = direct_match_details(record, "TP53")
    assert match and match["gene_matches"] == ["TP53"]


# 3
def test_direct_gene_synonym_match():
    record = base_record(gene="TP53")
    match = direct_match_details(record, "TP53S")
    assert match and "TP53S" in match["gene_matches"]


# 4
def test_direct_protein_name_match():
    record = base_record(name="Cellular tumor antigen p53")
    match = direct_match_details(record, "Cellular tumor antigen p53")
    assert match and match["protein_matches"] == ["Cellular tumor antigen p53"]


# 5
def test_discovery_does_not_combine_separate_fields():
    record = base_record()
    record["comments"] = [
        {"commentType": "FUNCTION", "texts": [{"value": "brain development"}]},
        {"commentType": "DISEASE", "texts": [{"value": "colorectal cancer"}]},
    ]
    assert annotation_match_contexts(record, "brain cancer") == []


# 6
def test_discovery_accepts_complete_phrase_in_one_field():
    record = base_record()
    record["comments"] = [{"commentType": "FUNCTION", "texts": [{"value": "Participates in DNA repair after damage."}]}]
    contexts = annotation_match_contexts(record, "DNA repair")
    assert len(contexts) == 1
    assert "DNA repair" in contexts[0]["text"]


# 7
def test_literature_title_alone_must_identify_record():
    record = base_record(gene="ABC1")
    record["references"] = [{"citation": {"title": "DNA repair in cancer"}}]
    assert annotation_match_contexts(record, "DNA repair") == []


# 8
def test_reviewed_status_is_swissprot():
    status = review_status(base_record(entry_type="UniProtKB reviewed (Swiss-Prot)"))
    assert status["status"] == "reviewed"
    assert "Swiss-Prot" in status["label"]


# 9
def test_unreviewed_status_is_trembl():
    status = review_status(base_record(entry_type="UniProtKB unreviewed (TrEMBL)"))
    assert status["status"] == "unreviewed"
    assert "TrEMBL" in status["label"]


# 10
def test_isoforms_are_read_from_alternative_products():
    record = base_record()
    record["comments"] = [{
        "commentType": "ALTERNATIVE PRODUCTS",
        "events": ["Alternative splicing"],
        "isoforms": [
            {"name": {"value": "1"}, "isoformIds": ["P1-1"], "isoformSequenceStatus": "Displayed"},
            {"name": {"value": "2"}, "isoformIds": ["P1-2"], "isoformSequenceStatus": "Described"},
        ],
    }]
    isoforms = extract_isoforms(record)
    assert [item["ids"][0] for item in isoforms] == ["P1-1", "P1-2"]


# 11
def test_isoform_explanation_separates_isoforms_from_accessions():
    assert "not the same thing as a separate UniProt accession" in ISOFORM_EXPLANATION


# 12
def test_experimental_evidence_code_is_labelled_conservatively():
    item = evidence_label({"evidenceCode": "ECO:0000269", "source": "PubMed", "id": "123"})
    assert item["category"] == "experimental"
    assert item["source_id"] == "123"


# 13
def test_similarity_evidence_code_is_not_called_experimental():
    item = evidence_label({"evidenceCode": "ECO:0000250", "source": "UniProtKB", "id": "Q1"})
    assert item["category"] == "similarity"
    assert "similarity" in item["label"].lower()


# 14
def test_function_evidence_reports_experimental_when_present():
    record = base_record()
    record["comments"] = [{
        "commentType": "FUNCTION",
        "texts": [{"value": "Does something", "evidences": [{"evidenceCode": "ECO:0000269"}]}],
    }]
    profile = function_evidence(record)
    assert "experimental evidence" in profile["summary"].lower()


# 15
def test_exact_sequence_duplicates_are_computed_only_from_sequences():
    one = base_record("P1")
    two = base_record("P2")
    three = base_record("P3")
    three["sequence"] = {"length": 100, "value": "B" * 100}
    duplicates = group_sequence_duplicates([one, two, three])
    assert duplicates["P1"] == ["P2"]
    assert "P3" not in duplicates


# 16
def test_domain_comparison_is_generic_not_vwfa_hardcoded():
    raw1 = base_record("P1")
    raw2 = base_record("P2")
    raw3 = base_record("P3")
    external = {"ensembl_canonical": None, "gene_centric_accession": None, "appris": {}}
    comparison = build_comparison_records([raw1, raw2, raw3], external)
    comparison[0]["features"]["domains"] = [{"name": "Kinase domain", "position": "10–80"}]
    comparison[1]["features"]["domains"] = [{"name": "Kinase domain", "position": "11–82"}]
    comparison[2]["features"]["domains"] = [{"name": "SH2 domain", "position": "20–70"}]
    explained = build_explanations(comparison)
    third_text = " ".join(explained[2]["sentences"])
    assert "Kinase domain" in third_text
    assert "VWFA" not in third_text


# 17
def test_no_domain_comparison_when_no_domain_is_common():
    raw1 = base_record("P1")
    raw2 = base_record("P2")
    external = {"ensembl_canonical": None, "gene_centric_accession": None, "appris": {}}
    comparison = build_comparison_records([raw1, raw2], external)
    comparison[0]["features"]["domains"] = [{"name": "Kinase domain", "position": None}]
    comparison[1]["features"]["domains"] = [{"name": "SH2 domain", "position": None}]
    explained = build_explanations(comparison)
    text = " ".join(sentence for item in explained for sentence in item["sentences"])
    assert "while it is annotated on" not in text


# 18
def test_result_links_include_uniprot_alphafold_and_pdb_search():
    links = external_links(base_record("P04637"))
    assert links["uniprot"].endswith("/P04637/entry")
    assert links["alphafold"].endswith("/P04637")
    assert links["pdb_search"].startswith("https://www.rcsb.org/search?request=")


# 19
def test_gene_record_retrieval_uses_stream_once_when_exact_succeeds():
    client = StreamClient(exact=[base_record("P1")])
    records, warnings = fetch_all_records_for_gene("TP53", 9606, client)
    assert len(records) == 1
    assert len(client.calls) == 1
    assert client.calls[0].startswith("(gene_exact:TP53)")
    assert warnings == []


# 20
def test_taxonomy_candidates_are_conservative_and_wordnet_free():
    candidates = candidate_species_phrases("TP53 in humans")
    assert "human" in candidates or "humans" in candidates
    # No semantic guess such as primate/mammal is manufactured locally.
    assert "primate" not in candidates
    assert "mammal" not in candidates


# Additional guard: placeholder taxa are not auto-selected.
def test_placeholder_taxonomy_record_is_detected():
    assert taxonomy_record_is_placeholder({"scientificName": "uncultured bacterium", "taxonId": 1}) is True

# 22 — plural organism wording is preserved for later removal from the query.
def test_species_resolution_preserves_plural_input_phrase():
    from app.core.taxonomy import resolve_species

    class TaxonomyClient:
        def concurrent_taxonomy_search(self, terms, size=25):
            human = {
                "scientificName": "Homo sapiens",
                "commonName": "Human",
                "taxonId": 9606,
                "rank": "species",
            }
            return {term: ([human] if term == "human" else []) for term in terms}, []

        def uniprot_search(self, query, size=30):
            return []

    result = resolve_species("anthrax in humans", TaxonomyClient())
    assert result["status"] == "resolved"
    assert result["match"]["record"]["taxonId"] == 9606
    assert result["match"]["input_phrase"] == "humans"


# 23 — the parser restores the Colab behavior: organism out, concept intact.
def test_natural_language_parser_extracts_anthrax_from_humans_query():
    from app.core.query import extract_search_concepts

    parsed = extract_search_concepts("anthrax in humans", "humans")
    assert parsed["primary"] == "anthrax"
    assert parsed["strict_variants"] == ["anthrax"]


# 24 — main lookup does not trap concept discovery inside the mentioned species.
def test_lookup_routes_anthrax_in_humans_to_global_discovery():
    from app.core.service import ProteinService

    anthrax_record = base_record(
        accession="A0AANTH",
        gene="pagA",
        name="Anthrax toxin protective antigen",
        entry_type="UniProtKB reviewed (Swiss-Prot)",
    )
    anthrax_record["organism"] = {
        "scientificName": "Bacillus anthracis",
        "commonName": None,
        "taxonId": 1392,
    }

    class SmartClient:
        def concurrent_taxonomy_search(self, terms, size=25):
            human = {
                "scientificName": "Homo sapiens",
                "commonName": "Human",
                "taxonId": 9606,
                "rank": "species",
            }
            return {term: ([human] if term == "human" else []) for term in terms}, []

        def concurrent_uniprot_search(self, queries, size=100):
            rows = {}
            for query in queries:
                rows[query] = [anthrax_record] if query == "(anthrax)" else []
            return rows, []

        def uniprot_search(self, query, size=30):
            return []

    result = ProteinService(SmartClient()).lookup("anthrax in humans")
    assert result["status"] == "no_direct_match"
    assert result["query"] == "anthrax"
    assert result["organism"]["taxon_id"] == 9606
    assert result["discovery"][0]["accession"] == "A0AANTH"
    assert result["discovery"][0]["bucket"] == "global"
    assert "Anthrax" in result["discovery"][0]["why"][0]["text"]


# 25 — discovery evidence may legitimately be a protein name, as in Colab v16.
def test_discovery_can_be_explained_by_protein_name():
    record = base_record(name="Anthrax toxin protective antigen")
    contexts = annotation_match_contexts(record, "anthrax")
    assert contexts
    assert contexts[0]["source"] == "Protein name"

# 26 — preserve Colab v16 direct-name behavior: a phrase inside a true protein name is direct identity evidence.
def test_partial_phrase_inside_protein_name_is_direct_identity():
    record = base_record(gene="ANTXR2", name="Anthrax toxin receptor 2")
    match = direct_match_details(record, "anthrax")
    assert match is not None
    assert match["gene"] == "ANTXR2"
    assert any("protein name" in reason.lower() for reason in match["reasons"])


# 27 — preserve the actual UniProt field name used by the Colab prototype.
def test_gene_labels_include_ordered_locus_names_plural_field():
    from app.core.records import get_all_gene_labels

    record = base_record()
    record["genes"][0]["orderedLocusNames"] = [{"value": "LOC123"}]
    assert "LOC123" in get_all_gene_labels(record)


# 28 — protein-name identity search includes short names and component names, as in Colab v16.
def test_protein_search_names_include_short_and_component_names():
    from app.core.records import get_all_protein_search_names

    record = base_record()
    record["proteinDescription"] = {
        "recommendedName": {
            "fullName": {"value": "Long protein name"},
            "shortNames": [{"value": "LPN"}],
        },
        "includes": [
            {"recommendedName": {"fullName": {"value": "Catalytic component"}}}
        ],
    }
    names = get_all_protein_search_names(record)
    assert "LPN" in names
    assert "Catalytic component" in names


# 29 — mentioned-organism discovery is grouped by gene and remains inspectable as the full gene record set.
def test_source_discovery_groups_by_gene_and_exposes_full_record_inspection():
    from app.core.search import discovery_search

    r1 = base_record("P1", gene="ANTXR2", name="Anthrax toxin receptor 2")
    r2 = base_record("P2", gene="ANTXR2", name="Anthrax toxin receptor 2")
    r1["comments"] = [{"commentType": "FUNCTION", "texts": [{"value": "Receptor involved in anthrax toxin uptake"}]}]
    r2["keywords"] = [{"name": "Anthrax"}]

    class Client:
        def concurrent_uniprot_search(self, queries, size=100):
            rows = {}
            for q in queries:
                if "organism_id:9606" in q and "virus_host_id" not in q:
                    rows[q] = [r1, r2]
                else:
                    rows[q] = []
            return rows, []

    payload = discovery_search("anthrax", Client(), taxon_id=9606, organism_name="Human")
    source = [item for item in payload["results"] if item["bucket"] == "mentioned_organism"]
    assert len(source) == 1
    assert source[0]["gene"] == "ANTXR2"
    assert source[0]["inspect"]["mode"] == "gene"
    assert source[0]["inspect"]["taxon_id"] == 9606
    assert source[0]["evidence_record_count"] == 2


# 30 — the main web front door restores the Colab behavior for "anthrax in humans":
# direct protein-name matches in Human are offered as proteins to inspect.
def test_lookup_anthrax_in_humans_returns_direct_protein_choices_when_name_fields_match():
    from app.core.service import ProteinService

    antxr1 = base_record("Q9H6X2", gene="ANTXR1", name="Anthrax toxin receptor 1")
    antxr2 = base_record("P58335", gene="ANTXR2", name="Anthrax toxin receptor 2")

    class Client:
        def concurrent_taxonomy_search(self, terms, size=25):
            human = {"scientificName": "Homo sapiens", "commonName": "Human", "taxonId": 9606, "rank": "species"}
            return {term: ([human] if term == "human" else []) for term in terms}, []

        def concurrent_uniprot_search(self, queries, size=100):
            rows = {}
            for q in queries:
                if "organism_id:9606" in q and "anthrax" in q:
                    rows[q] = [antxr1, antxr2]
                else:
                    rows[q] = []
            return rows, []

        def uniprot_search(self, query, size=30):
            return []

    result = ProteinService(Client()).lookup("anthrax in humans")
    assert result["status"] == "needs_protein_choice"
    assert {item["gene"] for item in result["options"]} == {"ANTXR1", "ANTXR2"}


# 31 — without an organism, the tool does not silently choose one species even if a name match exists globally.
def test_no_organism_does_not_silently_promote_global_direct_match():
    from app.core.service import ProteinService

    tp53 = base_record("P04637", gene="TP53", name="Cellular tumor antigen p53")

    class Client:
        def concurrent_taxonomy_search(self, terms, size=25):
            return {term: [] for term in terms}, []

        def concurrent_uniprot_search(self, queries, size=100):
            return {q: ([tp53] if "tp53" in q.lower() else []) for q in queries}, []

        def uniprot_search(self, query, size=30):
            return []

    result = ProteinService(Client()).lookup("TP53")
    assert result["status"] == "no_direct_match"
    assert result["organism"] is None
    assert result["discovery"]


# 32 — global discovery cards keep the Colab selection behavior: gene+source taxon opens all records.
def test_global_discovery_exposes_gene_level_inspection_when_possible():
    from app.core.search import discovery_search

    record = base_record("P04637", gene="TP53", name="Cellular tumor antigen p53")
    record["keywords"] = [{"name": "Li-Fraumeni syndrome"}]

    class Client:
        def concurrent_uniprot_search(self, queries, size=100):
            return {q: [record] for q in queries}, []

    payload = discovery_search("Li-Fraumeni syndrome", Client())
    item = payload["results"][0]
    assert item["bucket"] == "global"
    assert item["inspect"]["mode"] == "gene"
    assert item["inspect"]["gene"] == "TP53"
    assert item["inspect"]["taxon_id"] == 9606

# 33 — a direct identity result should not wait for global/virus discovery.
def test_direct_lookup_skips_unneeded_discovery_requests():
    from app.core.service import ProteinService

    tp53 = base_record("P04637", gene="TP53", name="Cellular tumor antigen p53")

    class Client:
        def __init__(self):
            self.search_calls = []
            self.stream_calls = []

        def concurrent_uniprot_search(self, queries, size=100):
            queries = list(queries)
            self.search_calls.extend(queries)
            return {q: ([tp53] if "organism_id:9606" in q else []) for q in queries}, []

        def uniprot_stream(self, query):
            self.stream_calls.append(query)
            return [tp53]

    client = Client()
    result = ProteinService(client).lookup(
        "TP53 in human",
        taxon_id=9606,
        organism_name="Human",
        organism_phrase="human",
        resolve_organism=False,
    )
    assert result["status"] == "ready"
    assert len(client.search_calls) == 1
    assert "organism_id:9606" in client.search_calls[0]
    assert not any("virus_host_id" in query for query in client.search_calls)
    assert not any(query == "(tp53)" for query in client.search_calls)
    assert len(client.stream_calls) == 1


# 34 — when direct identity fails, source rows are reused instead of searched twice.
def test_discovery_reuses_mentioned_organism_search_rows():
    from app.core.service import ProteinService

    receptor = base_record("P58335", gene="ANTXR2", name="Receptor protein")
    receptor["comments"] = [{"commentType": "FUNCTION", "texts": [{"value": "Receptor involved in anthrax toxin uptake"}]}]

    class Client:
        def __init__(self):
            self.calls = []

        def concurrent_uniprot_search(self, queries, size=100):
            queries = list(queries)
            self.calls.append(queries)
            out = {}
            for q in queries:
                if "organism_id:9606" in q and "virus_host_id" not in q:
                    out[q] = [receptor]
                else:
                    out[q] = []
            return out, []

    client = Client()
    result = ProteinService(client).lookup(
        "anthrax in humans",
        taxon_id=9606,
        organism_name="Human",
        organism_phrase="humans",
        resolve_organism=False,
    )
    assert result["status"] == "no_direct_match"
    flattened = [q for batch in client.calls for q in batch]
    source_queries = [q for q in flattened if "organism_id:9606" in q and "virus_host_id" not in q]
    assert len(source_queries) == 1


# 35 — UniProt explanations render without blocking on Ensembl/APPRIS by default.
def test_explain_records_defers_external_annotation_requests():
    from app.core.explain import explain_records

    class NoExternalClient:
        def gene_centric(self, *args, **kwargs):
            raise AssertionError("gene-centric should be lazy")
        def ensembl_lookup(self, *args, **kwargs):
            raise AssertionError("Ensembl should be lazy")
        def appris(self, *args, **kwargs):
            raise AssertionError("APPRIS should be lazy")

    payload = explain_records([base_record("P1")], "ABC1", "Human", NoExternalClient())
    assert payload["record_count"] == 1
    assert payload["external_annotations_loaded"] is False


# 36 — ready gene results expose enough metadata for bulk/selected downloads.
def test_explain_gene_exposes_download_scope_and_accessions():
    from app.core.service import ProteinService

    one = base_record("P1", gene="ABC1")
    two = base_record("P2", gene="ABC1", entry_type="UniProtKB unreviewed (TrEMBL)")

    class Client:
        def uniprot_stream(self, query):
            return [one, two]

    result = ProteinService(Client()).explain_gene("ABC1", 9606, species_name="Human")
    assert result["download_scope"] == {"mode": "gene", "gene": "ABC1", "taxon_id": 9606}
    assert result["accessions"] == ["P1", "P2"]
    assert result["reviewed_count"] == 1
    assert result["unreviewed_count"] == 1


# 37 — the web download surface includes the practical UniProt formats.
def test_download_formats_cover_fasta_tsv_json_xml_and_text():
    from app.main import DOWNLOAD_FORMATS

    assert {"fasta", "tsv", "json", "xml", "txt"}.issubset(DOWNLOAD_FORMATS)

# 38 — HTTP resilience: a transient 503 is retried, the successful response is cached,
# and the same request does not hit the upstream service again.
def test_http_client_retries_503_then_recovers_and_caches_success(tmp_path, monkeypatch):
    import json
    import threading
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

    from app.core.http import DataClient

    monkeypatch.delenv("VERCEL", raising=False)

    class Handler(BaseHTTPRequestHandler):
        calls = 0

        def do_GET(self):
            type(self).calls += 1

            if type(self).calls == 1:
                self.send_response(503)
                self.send_header("Content-Type", "application/json")
                self.end_headers()
                self.wfile.write(b'{"error":"temporary outage"}')
                return

            payload = json.dumps({"results": [{"primaryAccession": "PTEST"}]}).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)

        def log_message(self, format, *args):
            # Keep the test output quiet.
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        client = DataClient(
            cache_name=str(tmp_path / "retry_cache"),
            cache_days=1,
            timeout=3,
        )

        retry = client.session.get_adapter("http://").max_retries
        assert retry.status == 1
        assert retry.backoff_factor > 0
        assert 503 in retry.status_forcelist

        url = f"http://127.0.0.1:{server.server_port}/records"

        first = client._get_json("test source", url)
        assert first["results"][0]["primaryAccession"] == "PTEST"
        assert Handler.calls == 2  # initial 503 + successful retry

        second = client._get_json("test source", url)
        assert second == first
        assert Handler.calls == 2  # served from cache; no third upstream call
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


# 39 — when an exact protein identity exists, longer protein names that merely
# contain the same word must not compete with it as direct identities.
def test_exact_identity_beats_longer_name_mentions():
    from app.core.search import direct_search_phrases, group_direct_matches

    insulin = base_record("P01308", gene="INS", name="Insulin")
    atp5mk = base_record("Q9H2K0", gene="ATP5MK", name="Diabetes-associated protein in insulin-sensitive tissues")
    baiap2 = base_record("Q9UQB8", gene="BAIAP2", name="Insulin receptor substrate p53/p58")

    class Client:
        def concurrent_uniprot_search(self, queries, size=100):
            return {q: [insulin, atp5mk, baiap2] for q in queries}, []

    payload = direct_search_phrases(["insulin"], Client(), taxon_id=9606)
    groups = group_direct_matches(payload["matches"])
    assert [group["gene"] for group in groups] == ["INS"]


# 40 — if there is no exact identity, preserve the Colab contained-name
# fallback used by discovery-like protein names such as anthrax toxin receptors.
def test_contained_name_fallback_remains_when_no_exact_identity_exists():
    from app.core.search import direct_search_phrases, group_direct_matches

    antxr1 = base_record("Q9H6X2", gene="ANTXR1", name="Anthrax toxin receptor 1")
    antxr2 = base_record("P58335", gene="ANTXR2", name="Anthrax toxin receptor 2")

    class Client:
        def concurrent_uniprot_search(self, queries, size=100):
            return {q: [antxr1, antxr2] for q in queries}, []

    payload = direct_search_phrases(["anthrax"], Client(), taxon_id=9606)
    groups = group_direct_matches(payload["matches"])
    assert {group["gene"] for group in groups} == {"ANTXR1", "ANTXR2"}


# 41 — the first-screen summary should explain observable differences across a multi-record set.
def test_difference_summary_surfaces_record_set_variation():
    from app.core.explain import build_difference_summary

    rows = [
        {"length": 90, "name": "Protein A", "transcript": "ENST1", "existence": "protein level", "duplicates": [], "isoforms": []},
        {"length": 110, "name": "Protein A", "transcript": None, "existence": "transcript level", "duplicates": ["P3"], "isoforms": [{"ids": ["P2-1"]}]},
        {"length": 110, "name": "Protein B", "transcript": "ENST2", "existence": "transcript level", "duplicates": ["P2"], "isoforms": []},
    ]
    notes = build_difference_summary(rows)
    joined = " ".join(notes).lower()
    assert "90" in joined and "110" in joined
    assert "protein-name" in joined
    assert "transcript" in joined

# 42 — the browser UI must not explode a record card into hundreds of individual PDB links.
def test_frontend_uses_compact_pdb_structure_link():
    from pathlib import Path

    source = (Path(__file__).parents[1] / "app" / "static" / "app.js").read_text()
    assert "PDB structures (" in source
    assert "PDB ${escapeHtml(x.id)}" not in source


# 43 — the record page keeps the comparison explanation and practical downloads on the first screen.
def test_frontend_keeps_comparison_summary_and_full_set_downloads_visible():
    from pathlib import Path

    source = (Path(__file__).parents[1] / "app" / "static" / "app.js").read_text()
    assert "What differs across these records" in source
    assert "Download full record set" in source
    assert "What does reviewed vs unreviewed mean?" in source
