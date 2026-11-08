"""Phase 21 ADIM 3: frozen 2x2 tolerance/index diagnostic."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "phase21_reviewer1"
V1 = ROOT / "results/phase17/PHASE17_ADIM0_V1_REPRODUCTION.csv"
V2 = ROOT / "results/phase17/PHASE17_ADIM3_V2_CERT_PRIMARY.csv"
COMP5 = ["F", "S", "U", "G_claim", "R"]


def adjust(frame: pd.DataFrame, pcols: list[str]) -> pd.DataFrame:
    frame = frame.copy()
    frame["p_joint"] = frame[pcols].max(axis=1)
    frame["q_BH"] = np.nan
    frame["q_BY"] = np.nan
    for _, idx in frame.groupby("model").groups.items():
        p = frame.loc[idx, "p_joint"].to_numpy(float)
        frame.loc[idx, "q_BH"] = multipletests(p, alpha=.05, method="fdr_bh")[1]
        frame.loc[idx, "q_BY"] = multipletests(p, alpha=.05, method="fdr_by")[1]
    frame["certified_BH"] = frame.q_BH <= .05
    frame["certified_BY"] = frame.q_BY <= .05
    return frame


def main() -> None:
    OUT.mkdir(exist_ok=True)
    v1 = pd.read_csv(V1).sort_values(["model", "feature"]).reset_index(drop=True)
    v2 = pd.read_csv(V2).sort_values(["model", "feature"]).reset_index(drop=True)
    assert v1[["model", "feature"]].equals(v2[["model", "feature"]])

    rows: list[pd.DataFrame] = []
    # Absolute v1 / WI5 is the frozen official result.
    a5 = v1.copy()
    a5["WI"] = a5[COMP5].min(axis=1)
    a5["bottleneck"] = a5[COMP5].idxmin(axis=1)
    a5["point_adequate"] = a5.WI >= .8
    a5["tolerance"] = "absolute_v1"
    a5["index"] = "WI5"
    rows.append(a5)

    # Absolute v1 / WI6: N is the frozen scale-aware non-degeneracy score;
    # all original five scores remain exactly v1 absolute.
    a6 = v1.copy()
    a6["N"] = v2.N
    a6["p_N"] = v2.p_N
    comps6 = COMP5 + ["N"]
    a6["WI"] = a6[comps6].min(axis=1)
    a6["bottleneck"] = a6[comps6].idxmin(axis=1)
    a6["point_adequate"] = a6.WI >= .8
    # Frozen v1 p_joint already represents recursive-IUT over the five leaves.
    a6["p_base5"] = a6.p_joint
    a6 = adjust(a6, ["p_base5", "p_N"])
    a6["tolerance"] = "absolute_v1"
    a6["index"] = "WI6"
    rows.append(a6)

    # Scale-aware v2 / WI5, using the frozen primary leaf p-values.
    s5 = v2.copy()
    s5["WI"] = s5[COMP5].min(axis=1)
    s5["bottleneck"] = s5[COMP5].idxmin(axis=1)
    s5["point_adequate"] = s5.WI >= .8
    s5 = adjust(s5, [f"p_{c}" for c in COMP5])
    s5["tolerance"] = "scaled_v2_k1_rho1"
    s5["index"] = "WI5"
    rows.append(s5)

    # Scale-aware v2 / WI6 is the frozen Phase 17 primary result.
    s6 = v2.copy()
    s6["WI"] = s6.WI6
    s6["bottleneck"] = s6.bottleneck6
    s6["point_adequate"] = s6.point_adequate6
    s6["tolerance"] = "scaled_v2_k1_rho1"
    s6["index"] = "WI6"
    rows.append(s6)

    feature = pd.concat(rows, ignore_index=True, sort=False)
    keep = ["model", "feature", "tolerance", "index"] + comps6 + [
        "WI", "bottleneck", "point_adequate", "p_joint", "q_BH", "q_BY",
        "certified_BH", "certified_BY"]
    feature = feature.reindex(columns=keep)

    # Reproduction gates: both frozen diagonal cells must be exact.
    got_a5 = feature[(feature.tolerance == "absolute_v1") & (feature["index"] == "WI5")]
    np.testing.assert_array_equal(got_a5.certified_BH, v1.certified_BH)
    np.testing.assert_array_equal(got_a5.certified_BY, v1.certified_BY)
    np.testing.assert_allclose(got_a5.p_joint, v1.p_joint, atol=0, rtol=0)
    got_s6 = feature[(feature.tolerance == "scaled_v2_k1_rho1") & (feature["index"] == "WI6")]
    np.testing.assert_array_equal(got_s6.certified_BH, v2.certified_BH)
    np.testing.assert_array_equal(got_s6.certified_BY, v2.certified_BY)
    np.testing.assert_allclose(got_s6.p_joint, v2.p_joint, atol=0, rtol=0)

    summary = (feature.groupby(["model", "tolerance", "index"], sort=False)
               .agg(n_features=("feature", "size"), point_adequate=("point_adequate", "sum"),
                    BH=("certified_BH", "sum"), BY=("certified_BY", "sum"))
               .reset_index())
    bn = (feature.groupby(["model", "tolerance", "index", "bottleneck"])
          .size().rename("n").reset_index())
    summary = summary.merge(
        bn.groupby(["model", "tolerance", "index"])
        .apply(lambda z: json.dumps(dict(zip(z.bottleneck, z.n)), sort_keys=True), include_groups=False)
        .rename("bottleneck_distribution").reset_index(),
        on=["model", "tolerance", "index"], how="left")

    feature.to_csv(OUT / "PHASE21_ADIM3_FOURCELL_FEATURES.csv", index=False)
    summary.to_csv(OUT / "PHASE21_ADIM3_FOURCELL.csv", index=False)
    md = """# Faz 21 ADIM 3 — 2×2 tolerans/endeks tanısı

Bu tanı saklı Faz 17 tablolarından üretilmiştir; yeni model fit'i veya bootstrap yoktur. `absolute_v1 × WI5` ve `scaled_v2(k=1,rho=1) × WI6(kappa=0.25)` hücreleri dondurulmuş sonuçlarla bit-düzeyinde yeniden doğrulandı. WI6'nın N bileşeni ölçeğe bağlı tanımıyla korunmuş, “absolute” etiketi yalnız ilk beş bileşenin v1 mutlak toleranslarını ifade etmektedir.

## Özet

```text
""" + summary.to_string(index=False) + """
```

BH/BY ailesi her modelin 27 özelliğidir. Bu çıktı tanısaldır; dondurulmuş v1 kararlarını değiştirmez.
"""
    (OUT / "PHASE21_ADIM3_FOURCELL.md").write_text(md, encoding="utf-8")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
