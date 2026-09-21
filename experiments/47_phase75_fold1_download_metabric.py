#!/usr/bin/env python3
"""Acquire a separate targeted METABRIC mRNA cache for the Phase 7.5 fold 1 panel."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil

import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / "data/external/metabric"
RAW = ROOT / "data/external/metabric_phase75_fold1"
DATA = ROOT / "results" / "phase75_fold1" / "models"
BASE = "https://www.cbioportal.org/api"


def digest(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    records = []
    for name in ["genes.json", "patient_clinical.json", "sample_clinical.json", "samples.json"]:
        source, target = OLD / name, RAW / name
        if not target.exists():
            try:
                os.link(source, target)
                method = "hardlink to preserved Phase 3 raw response"
            except OSError:
                shutil.copy2(source, target)
                method = "copy of preserved Phase 3 raw response"
        else:
            method = "existing Phase 7.5 fold 1 cache"
        if digest(source) != digest(target):
            raise RuntimeError(f"Static raw cache differs: {name}")
        records.append({"file": name, "source": str(source), "method": method, "sha256": digest(target)})

    selected = pd.read_parquet(DATA / "X_train_selected.parquet").columns
    ensembl = [column[4:].split(".")[0] for column in selected if column.startswith("rna_")]
    mapping_path = RAW / "ensembl_symbol_mapping.json"
    if not mapping_path.exists():
        response = requests.post(
            "https://mygene.info/v3/query",
            data={"q": ",".join(ensembl), "scopes": "ensembl.gene", "fields": "symbol,entrezgene", "species": "human"},
            timeout=120,
        )
        response.raise_for_status()
        response.json()
        mapping_path.write_bytes(response.content)
    annotations = json.loads(mapping_path.read_bytes())
    records.append(
        {
            "file": mapping_path.name,
            "url": "https://mygene.info/v3/query",
            "method": "POST or verified cache",
            "requested_ensembl_count": len(ensembl),
            "sha256": digest(mapping_path),
        }
    )
    genes = json.loads((RAW / "genes.json").read_bytes())
    lookup = {row["hugoGeneSymbol"]: row["entrezGeneId"] for row in genes}
    matched = sorted({lookup[row["symbol"]] for row in annotations if row.get("symbol") in lookup})
    for start in range(0, len(matched), 25):
        path = RAW / f"mrna_batch_{start // 25:02d}.json"
        batch = matched[start : start + 25]
        if not path.exists():
            response = requests.post(
                BASE + "/molecular-profiles/brca_metabric_mrna/molecular-data/fetch",
                params={"projection": "DETAILED"},
                json={"sampleListId": "brca_metabric_mrna", "entrezGeneIds": batch},
                timeout=120,
            )
            response.raise_for_status()
            response.json()
            path.write_bytes(response.content)
        rows = json.loads(path.read_bytes())
        records.append(
            {
                "file": path.name,
                "url": BASE + "/molecular-profiles/brca_metabric_mrna/molecular-data/fetch",
                "requested_gene_count": len(batch),
                "returned_rows": len(rows),
                "sha256": digest(path),
            }
        )
        print(f"METABRIC phase75_fold1 batch {start // 25 + 1}/{(len(matched) + 24) // 25}: genes={len(batch)} rows={len(rows)}", flush=True)
    provenance = {
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Phase 7.5 fold 1 selected TCGA RNA panel; separate targeted cache",
        "requested_rna_features": len(ensembl),
        "matched_entrez_genes": len(matched),
        "records": records,
        "patient_or_sample_identifiers_exported": False,
    }
    (RAW / "download_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(json.dumps({"requested_rna_features": len(ensembl), "matched_entrez_genes": len(matched)}, indent=2))


if __name__ == "__main__":
    main()
