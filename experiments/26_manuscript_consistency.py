"""Verify the manuscript ledger and draft against frozen result artifacts.

This is a read-only diagnostic apart from its JSON and Markdown reports.  It
does not fit a model, call an explainer, or alter any scientific result.
"""

from __future__ import annotations

import csv
import json
import re
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "RESULTS_LEDGER.md"
DRAFT = ROOT / "MANUSCRIPT_METHODS_TO_AVAILABILITY.md"
OUT_DIR = ROOT / "results" / "manuscript_consistency"
JSON_OUT = OUT_DIR / "consistency.json"
MD_OUT = ROOT / "MANUSCRIPT_CONSISTENCY_REPORT.md"
BASELINES = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer" / "baseline_results.json"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def rows(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def criterion(path: Path, key: str) -> int:
    return int(next(r["count"] for r in rows(path) if r["criterion"] == key))


def main() -> int:
    started = time.perf_counter()
    mismatches: list[str] = []
    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(condition), "detail": detail})
        if not condition:
            mismatches.append(f"{name}: {detail}")

    ledger = LEDGER.read_text(encoding="utf-8")
    draft = DRAFT.read_text(encoding="utf-8")
    source_started = time.perf_counter()

    p3 = read_json(ROOT / "results/phase3_brca/metadata.json")
    p4 = read_json(ROOT / "results/phase4_rsf/metadata.json")
    rsf50 = read_json(ROOT / "results/phase56_rsf_n50/metadata.json")
    deep_pilot = read_json(ROOT / "results/phase6_deepsurv_pilot/metadata.json")
    deep_rev = read_json(ROOT / "results/phase61_deepsurv_explainer_revision/metadata.json")
    deep50 = read_json(ROOT / "results/phase62_deepsurv_n50/metadata.json")
    comparison = read_json(ROOT / "results/phase62_three_model_comparison/metadata.json")
    baselines = read_json(BASELINES)

    cox_summary = ROOT / "results/phase53_cox_recursive/summary.csv"
    rsf_summary = ROOT / "results/phase56_rsf_n50_certification/summary.csv"
    deep_summary = ROOT / "results/phase62_deepsurv_n50_certification/summary.csv"
    rsf_gate = rows(ROOT / "results/phase56_rsf_n50_null_gate/gate_summary.csv")[0]
    deep_gate = rows(ROOT / "results/phase62_deepsurv_n50_null_gate/gate_summary.csv")[0]
    margins = {r["model"]: r for r in rows(ROOT / "results/phase62_three_model_comparison/margin_distribution.csv")}
    r_delta = {r["scope"]: r for r in rows(ROOT / "results/phase56_comparison/R_delta_summary.csv")}
    distances = rows(ROOT / "results/phase55_rank_universe_diagnostic/distance_to_rsf.csv")

    values = {
        "train": p3["stats"]["original_train"],
        "test": p3["stats"]["original_test"],
        "selected": p3["stats"]["selected_features"],
        "cox_features": p3["stats"]["model_features"],
        "overlap": p3["stats"]["patient_overlap"],
        "audit": p3["stats"]["audit_unique_patients"],
        "external": p3["stats"]["external"]["complete_query_patients"],
        "panel": p3["stats"]["assessed_features"],
        "cox_n6": p3["stats"]["n6_point_adequate_features"],
        "cox_n50": p3["stats"]["point_adequate_features"],
        "cox_retained": p3["stats"]["previously_passing_still_passing"],
        "cox_bh": criterion(cox_summary, "recursive_studentized_BH"),
        "cox_by": criterion(cox_summary, "recursive_studentized_BY"),
        "rsf_fixed_n2": p4["results"]["point_adequate_n2_fixed_background_pilot"],
        "rsf_sampled_n2": p4["results"]["point_adequate_n2_sampled_background"],
        "rsf_n20": p4["results"]["point_adequate_n20_sampled_background"],
        "rsf_n50": rsf50["results"]["point_adequate_n50"],
        "rsf_bh": criterion(rsf_summary, "BH_certified_alpha_0_05"),
        "rsf_by": criterion(rsf_summary, "BY_certified_alpha_0_05"),
        "deep_n2": deep_rev["n2"]["point_adequate"],
        "deep_n50": deep50["results"]["point_adequate_n50"],
        "deep_bh": criterion(deep_summary, "BH_certified_alpha_0_05"),
        "deep_by": criterion(deep_summary, "BY_certified_alpha_0_05"),
    }

    expected = {
        "train": 673, "test": 169, "selected": 616, "cox_features": 614,
        "overlap": 23, "audit": 141, "external": 1466, "panel": 77,
        "cox_n6": 27, "cox_n50": 22, "cox_retained": 15, "cox_bh": 8, "cox_by": 8,
        "rsf_fixed_n2": 43, "rsf_sampled_n2": 35, "rsf_n20": 23, "rsf_n50": 20,
        "rsf_bh": 10, "rsf_by": 10, "deep_n2": 11, "deep_n50": 8,
        "deep_bh": 0, "deep_by": 0,
    }
    for name, value in values.items():
        check(f"artifact.{name}", value == expected[name], f"observed={value!r}, expected={expected[name]!r}")

    check("artifact.common_panel_composition", p3["stats"]["external"]["expression_genes"] == 72, "expected 72 RNA features")
    check("artifact.rsf_query_rows", rsf50["design"]["query_rows"] == {"TCGA": 12, "METABRIC": 12}, str(rsf50["design"]["query_rows"]))
    check("artifact.deep_query_rows", deep50["arguments"]["anchors_per_site"] == 12, str(deep50["arguments"]["anchors_per_site"]))
    check("artifact.cox_cindex", abs(baselines["models"]["coxph"]["test_c_index"] - 0.7064765300059418) < 1e-15, "baseline_results.json")
    check("artifact.rsf_cindex", abs(baselines["models"]["random_survival_forest"]["test_c_index"] - 0.714795008912656) < 1e-15, "baseline_results.json")
    check("artifact.deep_cindex", abs(deep_pilot["adim0"]["test_c_index"] - 0.7456922162804516) < 1e-15, "phase6 pilot metadata")
    check("artifact.rsf_gate", rsf_gate["gate_passed"] == "True" and rsf_gate["holm_flagged_controls"] == "0", str(rsf_gate))
    check("artifact.deep_gate", deep_gate["gate_passed"] == "True" and deep_gate["holm_flagged_controls"] == "0", str(deep_gate))
    check("artifact.deep_transition", comparison["transition_counts"] == {"neither": 59, "lost": 10, "new": 7, "retained": 1}, str(comparison["transition_counts"]))
    check("artifact.rsf_selected_delta", abs(float(r_delta["N20 point-pass 23"]["mean_delta_R"]) + 0.015689884649511904) < 1e-15, r_delta["N20 point-pass 23"]["mean_delta_R"])
    check("artifact.rank_universe_rejected", all(r["closer_after_restriction"] == "False" for r in distances), "all three distances must worsen")
    check("artifact.deep_margin", abs(float(margins["DeepSurv N=50"]["mean_margin"]) - 0.009305182599355491) < 1e-15, margins["DeepSurv N=50"]["mean_margin"])
    source_seconds = time.perf_counter() - source_started

    text_started = time.perf_counter()
    formula_block = """F_j = min_d (1/T) sum_t F_hat_{dtj}
  F_hat_{dtj} = 1 - (1/n) sum_i min(1, (phi_bar_{dtj}(x_i) - psi_tilde_j(f_dt,x_i))^2 / a_j^2)

S_j = min_d S_{dj}
  S_{dj} = 1 - (2/(T(T-1)n)) sum_{t<t'} sum_i min(1, (phi_bar_{dtj}(x_i)-phi_bar_{dt'j}(x_i))^2/a_j^2)

U_j = min_d U_{dj}
  U_{dj} = (1/T) sum_t [1 - min(1, (H_{dtj}-L_{dtj})/u_j)]

G_j = min_{d>=1} [1 - min(1, |I_hat_{dj} - I_hat_{0j}| / b_j)]

G_signed_j = min_d>=1 [1 - mean_i min(1, |phi_hat_dj(x_i)-phi_hat_0j(x_i)|/a_j)]
G_claim_j = min(G_j, G_signed_j)

R_j = min_d R_{dj}
  R_{dj} = 1 - (2/(T(T-1)(p-1))) sum_{t<t'} sum_{l!=j} |o_{dt,jl}-o_{dt',jl}|/2"""
    wi_block = """WI_j = min_k z_{jk}                              (Eq. 2.3.1)
Adequate_j = 1{z_jk > tau_k  for all k}           (Eq. 2.3.2)
WI_claim = min_{j in J*} WI_j                     (Eq. 2.3.3, J* önceden kayıtlı feature seti)"""
    iut = "Joint p-value: p_joint_j = max_k p_jk  (intersection-union test)"
    for name, block in (("component_formulas", formula_block), ("wi_formulas", wi_block), ("iut_formula", iut)):
        check(f"ledger.{name}", block in ledger, "verbatim block missing")
        check(f"draft.{name}", block in draft, "verbatim block missing")

    required_phrases = [
        "673-row training and 169-row test split with 616 selected inputs",
        "stored CoxPH model used 614 inputs",
        "yielding 141 disjoint TCGA audit patients",
        "contained 77 features: 72 RNA features",
        "1,466 METABRIC patients",
        "27/77 features passed at N=6, compared with 22/77 at N=50",
        "15 of the original 27 remained above threshold",
        "43/77 under fixed-mean background and 35/77 under sampled background",
        "23/77 features passed at N=20 and 20/77 at N=50",
        "11/77 features passed at N=2. At N=50, only 8/77 passed",
        "59 passed at neither depth, 10 were lost, seven were new at N=50, and only one was retained",
        "CoxPH N=50 | 22/77 | 8/77 | 8/77",
        "RSF N=50 | 20/77 | 10/77 | 10/77",
        "DeepSurv N=50 | 8/77 | 0/77 | 0/77",
        "0.250446",
        "0.001367",
        "0.00349",
        "0.743",
        "0.05857 to 0.06076",
        "0.02174 to 0.02381",
        "0/30 Holm-flagged violations",
        "4,999 bootstrap draws",
        "2,000 independent outer simulations, 999 inner resamples",
        "[TODO: Claude/kullanıcı ile netleştirilecek — public repo + lisans kararı bekliyor]",
    ]
    for phrase in required_phrases:
        check(f"draft.phrase.{phrase[:28]}", phrase in draft, f"missing: {phrase}")

    for forbidden in ("## Abstract", "## 1. Introduction", "## Related Work"):
        check(f"draft.forbidden.{forbidden}", forbidden not in draft, f"forbidden section present: {forbidden}")
    check("draft.no_external_citations", not re.search(r"\[[0-9]+\]|\([A-Z][A-Za-z-]+ et al\.,? \d{4}\)", draft), "citation-like token found")
    # Literal audit: every number printed in the manuscript must also occur in
    # the ledger, except subsection labels whose values are document structure.
    number_pattern = re.compile(
        r"(?<![\w.])[-+]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:[eE][-+]?\d+)?%?"
    )
    normalize_number = lambda value: value.replace(",", "").replace("−", "-")
    draft_numbers = {normalize_number(v) for v in number_pattern.findall(draft.replace("−", "-"))}
    ledger_numbers = {normalize_number(v) for v in number_pattern.findall(ledger.replace("−", "-"))}
    section_labels = {"2.1", "2.2", "2.3", "2.4", "3.1", "3.2", "3.3", "3.4", "4.1", "4.2", "4.3", "4.4"}
    uncovered_numbers = sorted(draft_numbers - ledger_numbers - section_labels)
    check(
        "draft.every_numeric_literal_in_ledger",
        not uncovered_numbers,
        f"uncovered numeric literals={uncovered_numbers}",
    )
    identifier_patterns = {
        "TCGA barcode": r"\bTCGA-[A-Z0-9]{2}-[A-Z0-9]{4}\b",
        "METABRIC sample ID": r"\bMB-[0-9]{4,}\b",
    }
    for label, pattern in identifier_patterns.items():
        hits = re.findall(pattern, ledger + "\n" + draft)
        check(f"privacy.{label}", not hits, f"found {len(hits)} identifier(s)")
    check("source.canonical_missing_disclosed", "section2-warrant-index-taslak.md` is absent" in ledger and "was unavailable" in draft, "canonical-source warning absent")
    text_seconds = time.perf_counter() - text_started

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    total_seconds = time.perf_counter() - started
    result = {
        "passed": not mismatches,
        "check_count": len(checks),
        "mismatch_count": len(mismatches),
        "mismatches": mismatches,
        "checks": checks,
        "timings_seconds": {
            "step1_source_and_ledger_validation": source_seconds,
            "step2_draft_validation": text_seconds,
            "step3_total_consistency_run": total_seconds,
        },
        "files_read": [
            str(LEDGER.relative_to(ROOT)), str(DRAFT.relative_to(ROOT)),
            "results/phase3_brca/metadata.json", "results/phase4_rsf/metadata.json",
            "results/phase56_rsf_n50/metadata.json", "results/phase6_deepsurv_pilot/metadata.json",
            "results/phase61_deepsurv_explainer_revision/metadata.json",
            "results/phase62_deepsurv_n50/metadata.json",
            "results/phase62_three_model_comparison/metadata.json",
            "results/phase53_cox_recursive/summary.csv",
            "results/phase56_rsf_n50_certification/summary.csv",
            "results/phase62_deepsurv_n50_certification/summary.csv",
            "results/phase56_rsf_n50_null_gate/gate_summary.csv",
            "results/phase62_deepsurv_n50_null_gate/gate_summary.csv",
            "results/phase56_comparison/R_delta_summary.csv",
            "results/phase55_rank_universe_diagnostic/distance_to_rsf.csv",
            "results/phase62_three_model_comparison/margin_distribution.csv",
            str(BASELINES),
        ],
        "files_written": [str(JSON_OUT.relative_to(ROOT)), str(MD_OUT.relative_to(ROOT))],
        "canonical_formula_source_present": False,
        "canonical_formula_source_note": "section2-warrant-index-taslak.md was not found; preserved Phase 2 specification was used and must be compared before submission.",
        "scientific_computation_performed": False,
    }
    JSON_OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    status = "PASS" if not mismatches else "FAIL"
    mismatch_md = "None." if not mismatches else "\n".join(f"- {m}" for m in mismatches)
    MD_OUT.write_text(
        f"""# Manuscript consistency report

Overall status: **{status}** ({len(checks)} checks; {len(mismatches)} mismatches).

## Step wall-clock

These are measured automated file-processing times. Editorial drafting time was not instrumented and is therefore not estimated.

| Step | Measured wall-clock (s) |
|---|---:|
| Step 1 — source re-read and ledger validation | {source_seconds:.6f} |
| Step 2 — manuscript/formula/scope/privacy validation | {text_seconds:.6f} |
| Step 3 — complete consistency program | {total_seconds:.6f} |

## Checks performed

- Re-read frozen JSON/CSV artifacts and the sibling baseline-results file.
- Compared all headline cohort, model, progression, gate, and certification counts used in the draft.
- Required byte-for-byte formula blocks in both ledger and draft.
- Checked prohibited sections, citation-like tokens, the availability placeholder, and common TCGA/METABRIC identifier patterns.
- Confirmed that this diagnostic performed no model or explainer calls.

## Mismatches

{mismatch_md}

## Unresolved source-provenance item

`section2-warrant-index-taslak.md` was not found. The ledger and draft therefore use the preserved Phase 2 specification plus the Phase 2.1 signed extension. This is disclosed in both files and must be checked against the canonical source before submission.

Machine-readable detail: `results/manuscript_consistency/consistency.json`.
""",
        encoding="utf-8",
    )
    print(json.dumps({"status": status, "checks": len(checks), "mismatches": len(mismatches), "seconds": total_seconds}))
    return 0 if not mismatches else 1


if __name__ == "__main__":
    raise SystemExit(main())
