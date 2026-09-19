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
