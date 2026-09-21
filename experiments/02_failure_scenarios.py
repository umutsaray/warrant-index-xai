"""Run eleven failure families plus a sign-flip claim regression (Section C.2)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from warrantxai.components import component_scores, generalizability_signed, importance, pairwise_orders
from warrantxai.scenarios import SCENARIOS, mechanism_shift
from _shared import CALIBRATION, markdown_table, save_run, score_rows


def run(seeds: int = 20, n: int = 1000, output: Path = Path("results")) -> pd.DataFrame:
    """Generate feature and scenario-level outcomes without tuning (Section C.3).

    Parameters
    ----------
    seeds, n : int
        Independent DGP repeats and shared-anchor sample size.
    output : Path
        Local artifact directory.

    Returns
    -------
    DataFrame
        All feature/component rows, including explicit unassessed G.
    """
    if seeds < 1:
        raise ValueError("seeds must be positive")
    started = perf_counter()
    rows = []
    checks = []
    notes = {}
    diagnostics = []
    for name, generator in SCENARIOS.items():
        for seed in range(seeds):
            s = generator(n=n, seed=20260916 + seed)
            z = component_scores(
                s.phi,
                s.reference,
                s.lower,
                s.upper,
                **{k: CALIBRATION[k] for k in ["a", "b", "u", "eta"]},
                common_support=s.common_support,
            )
            rows.extend(score_rows(z, scenario=name, seed=seed))
            notes[name] = s.assumptions
            imp = importance(s.phi)
            orders = pairwise_orders(imp, CALIBRATION["eta"])
            diagnostics.append(
                {
                    "scenario": name,
                    "seed": seed,
                    "non_tied_fraction": np.count_nonzero(orders)
                    / (np.prod(imp.shape[:-1]) * imp.shape[-1] * (imp.shape[-1] - 1)),
                    **s.diagnostics,
                }
            )
            if name == "oracle_invariant":
                passed = bool(np.min(z) > 0.8)
                expectation = "F/S/G/R=1; WI>0.8"
            elif name == "stable_but_wrong":
                passed = bool(np.mean(z[:, 0]) < 0.4 and np.min(z[:, [1, 4]]) > 0.99)
                expectation = "F mean<0.4, S/R>0.99"
            elif name == "noisy":
                passed = bool(np.mean(z[:, 2]) < 0.5)
                expectation = "U mean<0.5; S not penalized for isolated MC source"
            elif name == "training_instability":
                passed = bool(np.min(z[:, 1]) < 0.7 and np.min(z[:, 0]) == 1)
                expectation = "F=1 but S minimum<0.7"
            elif name == "covariate_shift":
                passed = bool(np.min(z[:, 3]) == 1)
                expectation = "standardized G=1"
            elif name == "mechanism_shift":
                passed = bool(z[0, 3] == 0 and np.min(z[:, 0]) == 1)
                expectation = "vanished feature G=0, F=1"
            elif name == "spurious_shortcut":
                passed = bool(z[-1, 3] < 0.2)
                expectation = "shortcut G<0.2"
            elif name == "near_tied_ranks":
                passed = bool(np.min(z[:, 4]) < 1)
                expectation = "some R<1; unique rank NOT certified"
            elif name == "redundant_duplicates":
                passed = bool(np.min(z[:2, 1]) < 0.9)
                expectation = "individual S<0.9, group sum invariant"
            elif name == "ood":
                passed = bool(np.isnan(z[:, 3]).all())
                expectation = "G/WI not assessable"
            else:
                passed = bool(np.all(z == 1))
                expectation = "zero model WI=1, not informative ranking"
            checks.append({"scenario": name, "seed": seed, "pattern_passed": passed, "expectation": expectation})
    frame = pd.DataFrame(rows)
    # mean skips NaN in pandas; preserve assessed-count and all-NA rows explicitly.
    summary = (
        frame.groupby("scenario", sort=False)[["F", "S", "U", "G", "R", "WI", "point_adequate"]].mean().reset_index()
    )
    summary["assessed_features"] = frame.groupby("scenario", sort=False).WI.count().to_numpy()
    check_frame = pd.DataFrame(checks)
    patterns = (
        check_frame.groupby(["scenario", "expectation"], sort=False).pattern_passed.agg(["sum", "count"]).reset_index()
    )
    flipped = mechanism_shift(n=n, seed=20260916, mode="sign_flip")
    zflip = component_scores(
        flipped.phi,
        flipped.reference,
        flipped.lower,
        flipped.upper,
        **{k: CALIBRATION[k] for k in ["a", "b", "u", "eta"]},
    )
    signed = generalizability_signed(flipped.phi.mean(axis=1), CALIBRATION["a"], flipped.common_support)
    magnitude = component_scores(
        flipped.phi,
        flipped.reference,
        flipped.lower,
        flipped.upper,
        **{k: CALIBRATION[k] for k in ["a", "b", "u", "eta"]},
        common_support=flipped.common_support,
        generalizability_mode="magnitude_only",
    )
    flip_rows = []
    for mode, scores in [("claim", zflip), ("magnitude_only", magnitude)]:
        mode_rows = score_rows(scores, scenario="mechanism_sign_flip", seed=20260916, generalizability_mode=mode)
        for j, row in enumerate(mode_rows):
            row.update(G_signed=signed[j], G_magnitude=magnitude[j, 3], G_claim=zflip[j, 3])
        flip_rows.extend(mode_rows)
    report = (
        "# Deney 02 — 11 failure scenario\n\n"
        f"Her senaryo: {seeds} bağımsız seed, n={n}, T=6, total sites=3, p=5. "
        "Tablo feature/seed ortalamalarıdır; bütün feature'ların adequate olduğu anlamına gelmez. "
        "G sütunu varsayılan G_claim değeridir. Point adequacy istatistiksel certification değildir.\n\n"
        + markdown_table(summary)
        + "\n\n## Önceden belirtilen davranış kontrolleri\n\n"
        + markdown_table(patterns)
        + "\n\n## Kör noktalar ve yorum\n\n"
        "- **Gerçek mechanism sign flip (Faz 2.1):** F/S/R=1; magnitude-only G=1. "
        f"Değişen feature G_signed={signed[0]:.6f}, G_claim={zflip[0, 3]:.6f}, WI={np.min(zflip[0]):.6f}. "
        f"WI>0.8: {bool(np.min(zflip[0]) > 0.8)}; magnitude-only sensitivity WI={np.min(magnitude[0]):.6f}. "
        "Formül değiştirilmeden verilen signed/claim tanımları uygulandı. "
        "assumption: phi_hat training-replicate ortalaması, clipping öncesinde alınır; küçük veya ortalamada iptal olan signed farklar hâlâ geçebilir.\n"
        "- **Near ties:** R ortalama pair agreement olduğu için belirsiz bir sıralama WI'yi eşik altına indirmeyebilir; unique rank garantisi yok.\n"
        "- **Zero model:** doğru sıfır açıklaması WI=1 alır; importance ya da model usefulness sonucu değildir.\n"
        "- **Noisy:** finite MC mean için bilinen noise modelinden ek CI genişliği kullanılır. Sonsuz MC expectation alınırsa hata kaybolur; bunu training instability diye göstermiyoruz.\n"
        "- **OOD:** NA değerleri düşük/başarılı puanlarla doldurulmadı; assessed coverage ayrıca gösterildi.\n\n"
        "## assumption kayıtları\n\n"
        + "\n".join(
            f"### {name}\n\n" + "\n".join(f"- {note}" for note in values) + "\n" for name, values in notes.items()
        )
    )
    save_run(output, "02_scenarios", frame, report, started, {"seeds": seeds, "n": n, "generalizability_mode": "claim"})
    check_frame.to_csv(output / "02_pattern_checks.csv", index=False)
    (output / "02_diagnostics.json").write_text(json.dumps(diagnostics, indent=2), encoding="utf-8")
    pd.DataFrame(flip_rows).to_csv(output / "02_sign_flip.csv", index=False)
    print(summary.to_string(index=False))
    print(patterns[["scenario", "sum", "count"]].to_string(index=False))
    return frame


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, default=20)
    parser.add_argument("--n", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    run(args.seeds, args.n, args.output)
