"""Phase 3 raw clinical and selected-panel mRNA REST download, no patient logs."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/external/metabric"
TCGA = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
BASE = "https://www.cbioportal.org/api"


def run() -> None:
    """Cache unmodified HTTP bodies and provenance (Phase 3 acquisition).

    Existing responses are reused verbatim. Targeted genes only, no bulk CNA,
    mutation or methylation profile. Sample keys never appear in logs.
    """
    RAW.mkdir(parents=True, exist_ok=True)
    records = []

    def get(
        name: str, url: str, *, params: dict | None = None, body: dict | None = None, form: dict | None = None
    ) -> list | dict:
        path = RAW / (name + ".json")
        existed = path.exists()
        if not existed:
            response = requests.request(
                "POST" if body is not None or form is not None else "GET",
                url,
                params=params,
                json=body,
                data=form,
                timeout=120,
            )
            if response.status_code != 200:
                raise RuntimeError(f"Request failed: resource={name}, HTTP={response.status_code}")
            response.json()  # validate before saving, never print response data
            path.write_bytes(response.content)
        data = path.read_bytes()
        records.append(
            {
                "file": path.name,
                "url": url,
                "params": params,
                "method": "POST" if body is not None or form is not None else "GET",
                "request": body or form,
                "bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(),
                "used_existing_raw_response": existed,
            }
        )
        result = json.loads(data)
        print(f"Raw resource {name}: records={len(result)}, bytes={len(data)}", flush=True)
        return result

    get("profiles", BASE + "/studies/brca_metabric/molecular-profiles")
    get("clinical_attributes", BASE + "/studies/brca_metabric/clinical-attributes")
    for kind in ["patient", "sample"]:
        get(
            kind + "_clinical",
            BASE + "/studies/brca_metabric/clinical-data",
            params={"clinicalDataType": kind.upper(), "projection": "DETAILED", "pageSize": 1000000},
        )
    get("samples", BASE + "/studies/brca_metabric/samples", params={"pageSize": 1000000})
    genes = get("genes", BASE + "/genes", params={"pageSize": 1000000, "projection": "DETAILED"})
    get("api_schema", BASE + "/v3/api-docs")
    get("sample_lists", BASE + "/studies/brca_metabric/sample-lists")
    cols = pd.read_parquet(TCGA / "X_train_selected.parquet").columns
    ids = [c[4:].split(".")[0] for c in cols if c.startswith("rna_")]
    annotations = get(
        "ensembl_symbol_mapping",
        "https://mygene.info/v3/query",
        form={"q": ",".join(ids), "scopes": "ensembl.gene", "fields": "symbol,entrezgene", "species": "human"},
    )
    lookup = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in genes}
    matched = sorted({lookup[a["symbol"]] for a in annotations if a.get("symbol") in lookup})
    for start in range(0, len(matched), 25):
        get(
            f"mrna_batch_{start // 25:02d}",
            BASE + "/molecular-profiles/brca_metabric_mrna/molecular-data/fetch",
            params={"projection": "DETAILED"},
            body={"sampleListId": "brca_metabric_mrna", "entrezGeneIds": matched[start : start + 25]},
        )
    provenance = {
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Unmodified raw HTTP response bodies; selected TCGA RNA panel only, not whole transcriptome",
        "records": records,
    }
    (RAW / "download_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")


if __name__ == "__main__":
    run()
