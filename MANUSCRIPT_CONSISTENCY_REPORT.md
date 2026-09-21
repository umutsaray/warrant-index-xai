# Manuscript consistency report

Overall status: **PASS** (73 checks; 0 mismatches).

## Step wall-clock

These are measured automated file-processing times. Editorial drafting time was not instrumented and is therefore not estimated.

| Step | Measured wall-clock (s) |
|---|---:|
| Step 1 — source re-read and ledger validation | 0.002039 |
| Step 2 — manuscript/formula/scope/privacy validation | 0.001835 |
| Step 3 — complete consistency program | 0.004168 |

## Checks performed

- Re-read frozen JSON/CSV artifacts and the sibling baseline-results file.
- Compared all headline cohort, model, progression, gate, and certification counts used in the draft.
- Required byte-for-byte formula blocks in both ledger and draft.
- Checked prohibited sections, citation-like tokens, the availability placeholder, and common TCGA/METABRIC identifier patterns.
- Confirmed that this diagnostic performed no model or explainer calls.

## Mismatches

None.

## Unresolved source-provenance item

`section2-warrant-index-taslak.md` was not found. The ledger and draft therefore use the preserved Phase 2 specification plus the Phase 2.1 signed extension. This is disclosed in both files and must be checked against the canonical source before submission.

Machine-readable detail: `results/manuscript_consistency/consistency.json`.
