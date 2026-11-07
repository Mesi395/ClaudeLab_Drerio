"""Split-half labellers for the latent-class comparison (Card C, step 1).

For each axis that the published annotation claims (dedifferentiated, proliferating), split the
protocol's marker list into two disjoint halves (alternate genes of the sorted list) and run
protocol v1 (`zf.protocol.label_axis`, all rules, all cells) on each half. Halves with < 2 genes
in the data are skipped (the protocol needs >= 2). The per-cell calls hold barcodes, so they are
written only to the git-ignored data/derived/ on the compute host.

Usage: python scripts/split_half_calls.py   ->  data/derived/splithalf_v1.tsv.gz
"""

import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_protocol import ROOT, axis_genes  # noqa: E402

from zf import io, protocol  # noqa: E402

AXES = ["dedifferentiated", "proliferating"]


def main() -> None:
    a = io.load()
    controls = yaml.safe_load((ROOT / "config" / "axes.yaml").read_text())["controls"]
    fib = next(c for c in controls if c["name"] == "ambient_fibroblast")
    kw = dict(ambient_genes=[g["symbol"] for g in fib["genes"]])
    pre = protocol.prepare(a)
    genes = axis_genes()
    cols = {}
    for axis in AXES:
        present = sorted(g for g in genes[axis] if g in a.var_names)
        for h, half in (("half1", present[0::2]), ("half2", present[1::2])):
            print(axis, h, half, flush=True)
            if len(half) < 2:
                print(f"  skipped: {len(half)} gene(s)", flush=True)
                continue
            out = protocol.label_axis(a, half, precomputed=pre, **kw)
            cols[f"{axis}.{h}.call"] = out["call"].to_numpy()
    df = pd.DataFrame(cols, index=a.obs_names)
    df.index.name = "cell"
    dest = ROOT / "data" / "derived" / "splithalf_v1.tsv.gz"
    df.to_csv(dest, sep="\t")
    print("wrote", dest, df.shape, {c: df[c].value_counts().to_dict() for c in df}, flush=True)


if __name__ == "__main__":
    main()
