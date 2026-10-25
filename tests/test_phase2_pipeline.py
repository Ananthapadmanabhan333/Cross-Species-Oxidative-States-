"""End-to-end Phase 2 on the fake APIs (SYNTHETIC fixtures; no network)."""

import copy
import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from extremonet.analysis.inputs import InputIntegrityError, verify_phase1_run
from extremonet.analysis.phase2 import run_phase2
from extremonet.data_collection.http_cache import sha256_file
from extremonet.pipeline import run_pipeline

from .conftest import REPO_ROOT, FakeResponse, FakeSession, string_api_handler

KEGG_GENES = {"P00001": "dra:DR_T001", "P00002": "dra:DR_T002", "Q00003": "dra:DR_T003", "Q00004": "dra:DR_T004",
              "P10001": "hsa:1001", "P10002": "hsa:1002", "P10003": "hsa:1003", "P10006": "hsa:1006"}
KO_LINKS = {
    "dra": "dra:DR_T001\tko:K90001\ndra:DR_T004\tko:K90004\ndra:DR_T005\tko:K90004\ndra:DR_T009\tko:K99999\n",
    "hsa": "hsa:1001\tko:K91001\nhsa:1006\tko:K90004\n",
}
CURATED = """species_key\tstring_id\tuniprot\tgene\tevidence_type\teffect_on_role\tassay\tfinding\tpmid\tdoi\tcuration_basis\tcurated_on
deinococcus_radiodurans\t243230.DR_T001\tP00001\ttstA\tperturbation\tsupports\tdeletion\tsynthetic: mutant sensitive\t1\t10.0/x\tsynthetic\t2026-01-01
deinococcus_radiodurans\t243230.DR_T001\tP00001\ttstA\tperturbation\tcontradicts\tdeletion\tsynthetic: mutant more resistant\t2\t10.0/y\tsynthetic\t2026-01-01
deinococcus_radiodurans\t243230.DR_T002\tP00002\ttstB\tperturbation\tsupports\tdeletion\tsynthetic: mutant sensitive\t3\t10.0/z\tsynthetic\t2026-01-01
"""
UNIPROT_EVIDENCE = {
    # accession: ec, panther, pfam
    "P00001": ("1.15.1.1", "PTHR11111;", "PF00081;"),       # Mn/Fe-like synthetic SOD
    "P00002": ("1.11.1.6", "PTHR22222;", "PF00199;"),
    "Q00003": ("", "PTHR33333;", "PF00578;"),
    "Q00004": ("1.8.4.12", "PTHR44444;", "PF01641;"),
    "P10001": ("1.15.1.1", "PTHR99999;", "PF00080;"),       # Cu/Zn-like: same EC, no shared family
    "P10002": ("1.11.1.6", "PTHR22222;", "PF00199;"),       # homolog of P00002
    "P10003": ("1.11.1.24", "PTHR55555;", "PF08534;"),
    "P10005": ("", "PTHR66666;", "PF09999;"),
    "P10006": ("1.8.4.12", "PTHR77777;", "PF07777;"),
}
ANNOTATIONS = {  # string id -> list of (category, term)
    "243230.DR_T001": [("Process", "GO:0006801"), ("Function", "GO:0004784"), ("KEGG", "dra04146")],
    "243230.DR_T002": [("Process", "GO:0042744"), ("Function", "GO:0004096"), ("KEGG", "dra00380")],
    "243230.DR_T003": [("Process", "GO:0006979"), ("InterPro", "IPR036249")],
    "243230.DR_T004": [("Process", "GO:0006979")],
    "9606.ENSPTEST01": [("Process", "GO:0006801"), ("Function", "GO:0004784"), ("KEGG", "hsa04146")],
    "9606.ENSPTEST02": [("Process", "GO:0042744"), ("Function", "GO:0004096"), ("KEGG", "hsa00380")],
    "9606.ENSPTEST03": [("Process", "GO:0006979"), ("InterPro", "IPR036249")],
}


def phase2_handler(method, url, params):
    if url.endswith("/json/functional_annotation"):
        ids = params["identifiers"].split("\r")
        terms = {}
        for sid in ids:
            for cat, term in ANNOTATIONS.get(sid, []):
                terms.setdefault((cat, term), []).append(sid)
        rows = [{"category": c, "term": t, "inputGenes": g, "description": f"synthetic {t}"} for (c, t), g in terms.items()]
        return FakeResponse(body=json.dumps(rows))
    if url.endswith("/json/homology_best"):
        hits = {"243230.DR_T002": ("9606.ENSPTEST02", 300.0), "9606.ENSPTEST02": ("243230.DR_T002", 290.0)}
        rows = [{"ncbiTaxonId_A": 0, "stringId_A": q, "ncbiTaxonId_B": 0, "stringId_B": hits[q][0], "bitscore": hits[q][1]}
                for q in params["identifiers"].split("\r") if q in hits]
        return FakeResponse(body=json.dumps(rows))
    if "rest.kegg.jp/link/ko/" in url:
        return FakeResponse(body=KO_LINKS[url.rsplit("/", 1)[1]])
    if url.endswith("/json/get_string_ids") and params["identifiers"].split("\r")[0].startswith("DR_"):
        rows = [{"queryIndex": i, "queryItem": loc, "stringId": f"243230.{loc}", "ncbiTaxonId": 243230,
                 "taxonName": "synthetic", "preferredName": loc, "annotation": "synthetic"}
                for i, loc in enumerate(params["identifiers"].split("\r")) if loc != "DR_T009"]
        return FakeResponse(body=json.dumps(rows))
    if url.endswith("/json/network") and "243230.DR_T005" in params["identifiers"]:
        rows = json.loads(string_api_handler(method, url, params).text)
        rows.append({"stringId_A": "243230.DR_T004", "stringId_B": "243230.DR_T005", "preferredName_A": "a",
                     "preferredName_B": "b", "ncbiTaxonId": 243230, "score": 0.85, "nscore": 0, "fscore": 0,
                     "pscore": 0, "ascore": 0, "escore": 0, "dscore": 0, "tscore": 0.5})
        return FakeResponse(body=json.dumps(rows))
    if url.endswith("/json/ppi_enrichment"):
        return FakeResponse(body=json.dumps([{"number_of_nodes": 4, "number_of_edges": 2, "average_node_degree": 1,
                                              "local_clustering_coefficient": 0, "expected_number_of_edges": 0,
                                              "p_value": 0.01}]))
    if url.endswith("/uniprotkb/stream") and "ec" in params.get("fields", ""):
        lines = ["Entry\tEC number\tPANTHER\tPfam\tKEGG\tReactome\tPubMed ID"]
        for acc, (ec, panther, pfam) in UNIPROT_EVIDENCE.items():
            if ("1299" in params["query"]) == acc.startswith(("P0000", "Q0000")):
                lines.append("\t".join([acc, ec, panther, pfam, KEGG_GENES.get(acc, ""), "", "123"]))
        return FakeResponse(body="\n".join(lines) + "\n")
    if "QuickGO" in url:
        dr = params["taxonId"] == 1299
        results = ([{"geneProductId": "UniProtKB:P00001", "goId": "GO:0006801", "qualifier": "involved_in",
                     "goEvidence": "IDA", "evidenceCode": "ECO:0000314", "reference": "PMID:111"},
                    {"geneProductId": "UniProtKB:P00002", "goId": "GO:0042744", "qualifier": "involved_in",
                     "goEvidence": "IEA", "evidenceCode": "ECO:0000501", "reference": "GO_REF:1"}] if dr else
                   [{"geneProductId": "UniProtKB:P10001", "goId": "GO:0006801", "qualifier": "involved_in",
                     "goEvidence": "IMP", "evidenceCode": "ECO:0000315", "reference": "PMID:222"},
                    {"geneProductId": "UniProtKB:P10003", "goId": "GO:0006979", "qualifier": "NOT|involved_in",
                     "goEvidence": "IDA", "evidenceCode": "ECO:0000314", "reference": "PMID:333"}])
        return FakeResponse(body=json.dumps({"numberOfHits": len(results), "results": results,
                                             "pageInfo": {"resultsPerPage": 200, "current": 1, "total": 1}}))
    return string_api_handler(method, url, params)


@pytest.fixture
def phase2_setup(config, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    session = FakeSession(phase2_handler)
    p1 = run_pipeline(config, session=session, sleep=lambda s: None)
    cfg = yaml.safe_load((REPO_ROOT / "config" / "phase2_analysis.yaml").read_text())
    cfg["phase1"]["run_id"] = p1["run_id"]
    cfg["paths"] = {"results_dir": str(tmp_path / "results" / "phase2"), "manifests_dir": config["paths"]["manifests_dir"]}
    cfg["network"]["sensitivity_thresholds"] = [0.7, 0.9]
    cfg["communities"].update({"n_seeds": 3, "resolutions_for_stability": [1.0]})
    cfg["topology"]["size_matched_samples"] = 20
    cfg["alignment"].update({"null_permutations": 20, "rewired_null_samples": 3, "max_iterations": 50})
    cfg["ranking"]["sensitivity"]["dirichlet_samples"] = 10
    cfg["ranking"]["min_evidence_coverage"] = 0.2
    cfg["embedding"]["min_edges_for_evaluation"] = 1000
    cfg["coverage_variant"].update({"rewired_null_samples": 2, "max_iterations": 30})
    curated = tmp_path / "curated.tsv"
    curated.write_text(CURATED)
    cfg["curated_evidence"] = str(curated)
    cfg["_phase1"] = config
    return cfg, session, p1


def test_phase2_end_to_end(phase2_setup, tmp_path):
    cfg, session, p1 = phase2_setup
    result = run_phase2(copy.deepcopy(cfg), session=session, sleep=lambda s: None)
    out = result["results_dir"]
    for rel in ("technical_report.md", "analysis_summary.json", "manifest.json", "SHA256SUMS.txt",
                "candidates/candidate_pairs_all.tsv", "candidates/protein_evidence.tsv",
                "cross_species/alignment_summary.tsv", "cross_species/annotation_overlap_by_category.tsv",
                "evaluation/threshold_sensitivity_networks.tsv", "evaluation/similarity_measure_sensitivity.tsv",
                "deinococcus_radiodurans/node_metrics.tsv", "homo_sapiens/module_enrichment.tsv",
                "figures/degree_ccdf.png", "figures/ranking_weight_sensitivity.png"):
        assert (out / rel).is_file(), rel

    ranked = pd.read_csv(out / "candidates" / "candidate_pairs_all.tsv", sep="\t")
    assert len(ranked) == 4 * 5                                  # every D.rad x human node pair
    pair = ranked.set_index(["dr_string_id", "hs_string_id"])
    assert pair.loc[("243230.DR_T002", "9606.ENSPTEST02"), "category"] == "ortholog_rbh"
    # same EC, different families, shared KEGG map + experimental evidence -> candidate analogue
    assert pair.loc[("243230.DR_T001", "9606.ENSPTEST01"), "category"] == "candidate_analogue"
    # shared InterPro entry counts as (weak) homology, never as an analogue
    assert pair.loc[("243230.DR_T003", "9606.ENSPTEST03"), "homology_tier"] == "shared_interpro_entry"
    # NOT-qualified annotation is recorded as contradictory evidence
    assert "NOT-qualified" in pair.loc[("243230.DR_T003", "9606.ENSPTEST03"), "contradictory_evidence"]
    for col in ("supporting_evidence", "uncertainty", "proposed_validation", "references"):
        assert col in ranked
    # shared KEGG KO without a reciprocal best hit -> KO ortholog
    assert pair.loc[("243230.DR_T004", "9606.ENSPTEST06"), "category"] == "ortholog_ko"
    # curated mixed evidence (0.5) and GO IMP for the human protein (1.0): mean over available
    mixed = pair.loc[("243230.DR_T001", "9606.ENSPTEST01")]
    assert mixed["perturbation"] == pytest.approx(0.75)
    assert "curated perturbation evidence against" in mixed["contradictory_evidence"]
    # evidence for one protein only is not halved: missing evidence is not negative
    assert pair.loc[("243230.DR_T002", "9606.ENSPTEST03"), "perturbation"] == pytest.approx(1.0)
    assert pd.isna(pair.loc[("243230.DR_T004", "9606.ENSPTEST06"), "perturbation"])

    summary = json.loads((out / "analysis_summary.json").read_text())
    assert summary["embedding"]["deinococcus_radiodurans"]["evaluated"] is False
    assert summary["proxy_evaluation"]["components"]["ec"]["evaluable"] is False   # 2 positives < 10
    cov = summary["coverage_variant"]
    assert cov["enabled"] and cov["projected_dr_genes"] == 2 and cov["mapped_to_string"] == 5
    variant_nodes = pd.read_csv(out / "coverage_variant" / "nodes.tsv", sep="\t")
    assert set(variant_nodes.loc[~variant_nodes["is_original_seed"], "string_id"]) == {"243230.DR_T005"}
    assert summary["curated_evidence"]["rows"] == 3
    manifest = json.loads(result["manifest"].read_text())
    for item in manifest["outputs"]:
        assert sha256_file(tmp_path / item["path"]) == item["sha256"]
    assert {r["source"] for r in manifest["raw_inputs"]} == {"string", "uniprot", "quickgo", "kegg"}


def test_phase2_offline_rerun_is_identical(phase2_setup, offline_session):
    cfg, session, _ = phase2_setup
    first = run_phase2(copy.deepcopy(cfg), session=session, sleep=lambda s: None)
    offline = copy.deepcopy(cfg)
    offline["mode"] = "offline"
    offline["_phase1"]["mode"] = "offline"
    second = run_phase2(offline, session=offline_session, sleep=lambda s: None)
    assert offline_session.calls == []
    for rel in ("candidates/candidate_pairs_all.tsv", "deinococcus_radiodurans/node_metrics.tsv",
                "homo_sapiens/module_enrichment.tsv", "cross_species/alignment_summary.tsv"):
        a = pd.read_csv(first["results_dir"] / rel, sep="\t")
        b = pd.read_csv(second["results_dir"] / rel, sep="\t")
        pd.testing.assert_frame_equal(a, b)


def test_curated_evidence_must_reference_network_proteins(phase2_setup, tmp_path):
    cfg, session, _ = phase2_setup
    bad = tmp_path / "bad_curated.tsv"
    bad.write_text(CURATED.replace("243230.DR_T002", "243230.DR_NOPE"))
    cfg = copy.deepcopy(cfg)
    cfg["curated_evidence"] = str(bad)
    with pytest.raises(ValueError, match="DR_NOPE"):
        run_phase2(cfg, session=session, sleep=lambda s: None)


def test_phase2_refuses_tampered_phase1_inputs(phase2_setup, tmp_path):
    cfg, _, p1 = phase2_setup
    nodes = Path(p1["processed_dir"]) / "homo_sapiens" / "nodes.tsv"
    nodes.write_text(nodes.read_text() + "\n")
    with pytest.raises(InputIntegrityError):
        verify_phase1_run(Path(p1["manifest"]), tmp_path)
