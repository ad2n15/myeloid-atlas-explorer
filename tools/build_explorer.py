#!/usr/bin/env python3
"""Data files for the static atlas explorer (explorer/; GitHub Pages).

Per atlas (BM25, PBMC4) a stratified subsample of cells (level2 x study, rare
strata kept whole) is exported as compact little-endian binaries plus a
meta.json describing them (format: explorer/DATA_FORMAT.md):

  data/<atlas>/meta.json           fields, categories, colours, gene index, sizes
  data/<atlas>/umap.bin            global scVI UMAP, uint16 x,y interleaved
  data/<atlas>/l2umap/<lin>.bin    Level 2 lineage UMAP: uint32 row ids, then uint16 x,y
  data/<atlas>/obs/<field>.bin     categorical codes (uint8/uint16) or quantised numbers (uint16)
  data/<atlas>/genes/<gene>.bin.gz nonzero bitmap + uint8 log1p values (gzip)
  data/<atlas>/scores/<name>.bin.gz per-cell cytokine-response score, int8 z (-127 = not scored)
  data/results/*.json              donor-level results shared by both atlases

Per-cell scores are descriptive: for each cell whose Level 2 maps to a huCIRA
dictionary cell type, the signature score is the mean log1p of the signature
genes minus that of expression-matched control genes (AddModuleScore-style,
25 expression bins, 50 controls per gene), z-scored within the cell type.
Inference in the paper is donor-level.

  python3 scripts/build_explorer.py --atlas PBMC4 --n-cells 60000
  python3 scripts/build_explorer.py --atlas BM25 --n-cells 150000
  python3 scripts/build_explorer.py --results-only
  python3 scripts/build_explorer.py --donor-fields BM25    # refresh donor fields after a donor_map.csv change

Donor-level result tables carry `tested` (scripts/test_eligibility.py): a contrast x cell type is
reported and enters BH only if the donor pseudobulks hold a median of >= min_median_cells_tested
cells in each arm; untested rows keep their effect with q = NaN.
"""

from __future__ import annotations

import argparse
import colorsys
import gzip
import json
import shutil
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from annotate_level1 import atlas_h5ad  # noqa: E402
from pseudobulk_aggregate import cell_unit_keys, map_unit_keys  # noqa: E402
from response_modules import GROUP_LABEL, module_order, pretty, read_groups  # noqa: E402

OUT = REPO / "explorer/site/data"
GMT = REPO / "metadata/hucira_signatures/hcd_up_p0.01_lfc1.gmt"
SEED = 0
N_HVG = 2000
MIN_DETECT = 0.001
BUDGET_MB = {"BM25": 200, "PBMC4": 90}

# obs fields exported per atlas: (field id, label, source column, kind)
COMMON = [
    ("gse", "Study", "gse", "cat"),
    ("sample_id", "Sample", "sample_id", "cat"),
    ("gsm", "GEO sample (GSM)", "gsm", "cat"),
    ("disease_groups", "Disease group (study metadata)", "disease_groups", "cat"),
    ("primary_disease", "Primary disease (study metadata)", "primary_disease", "cat"),
    ("ct_bm_pred", "CellTypist bone marrow (predicted)", "predicted_labels_bone_marrow", "cat"),
    ("ct_bm_mv", "CellTypist bone marrow (majority voting)", "majority_voting_labels_bone_marrow", "cat"),
    ("ct_blood_pred", "CellTypist blood (predicted)", "predicted_labels_blood", "cat"),
    ("ct_blood_mv", "CellTypist blood (majority voting)", "majority_voting_labels_blood", "cat"),
    ("ct_low_pred", "CellTypist Immune_All_Low (predicted)", "predicted_labels_immune_low", "cat"),
    ("ct_low_mv", "CellTypist Immune_All_Low (majority voting)", "majority_voting_labels_immune_low", "cat"),
    ("ct_high_pred", "CellTypist Immune_All_High (predicted)", "predicted_labels_immune_high", "cat"),
    ("ct_high_mv", "CellTypist Immune_All_High (majority voting)", "majority_voting_labels_immune_high", "cat"),
    ("n_genes", "Genes detected", "n_genes_by_counts", "num"),
    ("total_counts", "UMI counts", "total_counts", "num"),
    ("pct_mt", "Mitochondrial reads (%)", "pct_counts_mt", "num"),
    ("doublet_score", "Scrublet doublet score", "doublet_score", "num"),
    ("ct_bm_conf", "CellTypist bone marrow confidence", "conf_score_bone_marrow", "num"),
]
EXTRA = {
    "BM25": [("author_state", "Author label (cell state)", "cell state", "cat"),
             ("author_annot", "Author label (GEO annotation)", "geo_Annotation", "cat"),
             ("age", "Age (years, where deposited)", "age", "num")],
    "PBMC4": [("author_type", "Author label (cell type)", "cell_type", "cat"),
              ("sctype", "Author label (scType)", "scType_celltype", "cat"),
              ("mutation_group", "CH mutation group", "MUTATION.GROUP", "cat"),
              ("clone", "Clone call (mutant / wild-type)", "Clone", "cat")],
}
DONOR_FIELDS = [("donor", "Donor"), ("condition", "Condition"), ("stage", "Stage"), ("timepoint", "Timepoint"),
                ("sex", "Sex"), ("genotype", "CH genotype"), ("tissue", "Tissue"), ("population", "Sorted population")]

# Level 1 lineage hues (fixed, so a lineage keeps its colour in both atlases)
LINEAGE_HUE = {
    "HSPC": "#2a78d6", "MkEry_prog": "#1baf7a", "Erythroid": "#c0392b", "Megakaryocyte": "#8e44ad",
    "Neutrophil_lineage": "#eda100", "Monocyte": "#eb6834", "cDC": "#e87ba4", "pDC": "#b5527f",
    "Mast": "#8c5a00", "B_prog": "#00a5b5", "B": "#0b7a75", "Plasma": "#5d6d7e", "T": "#6250d6",
    "NK": "#3f9c35", "Stromal": "#7f6a4d", "Endothelial": "#a04000", "Blast_like": "#d4ac0d",
    "Unresolved": "#a6a6a6", "Doublet": "#cccccc",
}
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#6250d6", "#e87ba4", "#00a5b5", "#8c5a00",
               "#c0392b", "#3f9c35", "#b5527f", "#5d6d7e"]


# ── helpers ──────────────────────────────────────────────────────────────────
def read_modules() -> pd.DataFrame:
    """Response modules (docs/phase1/response_modules.tsv) with their display rank (consensus-group order,
    response_modules.module_order) and whether the module is one of the huCIRA highly similar groups."""
    mods = pd.read_csv(REPO / "docs/phase1/response_modules.tsv", sep="\t")
    rank = {m: i for i, m in enumerate(module_order(mods))}
    return mods.assign(module_rank=mods["module"].map(rank),
                       module_grouped=mods["module"].isin(set(GROUP_LABEL.values())))


def module_of(cytokines: pd.Series) -> pd.Series:
    """Module of a cytokine by the fixed rule (same in every cell type): its highly similar group, else itself."""
    hs = read_groups()["highly_similar_group"].reindex(cytokines.values)
    return pd.Series([GROUP_LABEL[int(g)] if pd.notna(g) else pretty(c) for c, g in zip(cytokines, hs)],
                     index=cytokines.index)


def read_gmt(path: Path) -> dict[str, list[str]]:
    out = {}
    for line in path.read_text().splitlines():
        name, _, *genes = line.split("\t")
        out[name] = genes
    return out


def hex_shades(base: str, n: int) -> list[str]:
    """n distinguishable shades of one hue (lightness/saturation spread)."""
    r, g, b = (int(base[i:i + 2], 16) / 255 for i in (1, 3, 5))
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    if n == 1:
        return [base]
    out = []
    for k in range(n):
        lk = 0.32 + 0.38 * k / (n - 1)
        hk = (h + (k % 2) * 0.03) % 1
        rr, gg, bb = colorsys.hls_to_rgb(hk, lk, min(1, s * (1.05 - 0.25 * k / (n - 1))))
        out.append("#%02x%02x%02x" % (round(rr * 255), round(gg * 255), round(bb * 255)))
    return out


def palette(cats: list[str]) -> list[str]:
    """Fixed-order categorical palette; beyond 12 categories, golden-angle hues (marked in meta)."""
    if len(cats) <= len(CATEGORICAL):
        return CATEGORICAL[:len(cats)]
    out = []
    for k in range(len(cats)):
        h = (0.61 + k * 0.381966) % 1
        rr, gg, bb = colorsys.hls_to_rgb(h, 0.45 + 0.12 * (k % 3) / 2, 0.65)
        out.append("#%02x%02x%02x" % (round(rr * 255), round(gg * 255), round(bb * 255)))
    return out


def read_obs(f: h5py.File, col: str) -> pd.Series:
    return pd.Series(ad.io.read_elem(f["obs"][col]))


def quantise(v: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """float -> uint16 (0 = missing, 1..65535 = lo..hi)."""
    q = np.zeros(len(v), np.uint16)
    ok = np.isfinite(v)
    span = hi - lo if hi > lo else 1.0
    q[ok] = 1 + np.round((np.clip(v[ok], lo, hi) - lo) / span * 65534).astype(np.uint16)
    return q


def write(path: Path, arr: np.ndarray, gz: bool = False) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    b = np.ascontiguousarray(arr).astype(arr.dtype.newbyteorder("<"), copy=False).tobytes()
    if gz:
        b = gzip.compress(b, compresslevel=9, mtime=0)
    path.write_bytes(b)
    return len(b)


def subsample(strata: pd.Series, n: int, rng: np.random.Generator) -> np.ndarray:
    """Indices of a stratified subsample with a common per-stratum cap (small strata kept whole)."""
    sizes = strata.value_counts()
    lo, hi = 1, int(sizes.max())
    while lo < hi:                                     # largest cap with total <= n
        mid = (lo + hi + 1) // 2
        lo, hi = (mid, hi) if np.minimum(sizes, mid).sum() <= n else (lo, mid - 1)
    keep = []
    codes, uniq = pd.factorize(strata)
    for k in range(len(uniq)):
        idx = np.flatnonzero(codes == k)
        keep.append(idx if len(idx) <= lo else rng.choice(idx, lo, replace=False))
    return np.sort(np.concatenate(keep))


def read_rows(f: h5py.File, rows: np.ndarray, block: int = 50_000) -> sparse.csr_matrix:
    """Selected rows of the CSR matrix X (log1p) as float32 CSR, read block by block."""
    X = f["X"]
    n_obs, n_var = X.attrs["shape"]
    indptr = X["indptr"]
    datas, inds, lens = [], [], []
    for b0 in range(0, n_obs, block):
        sel = rows[(rows >= b0) & (rows < b0 + block)]
        if not len(sel):
            continue
        ptr = indptr[sel[0]:sel[-1] + 2].astype(np.int64)
        a = ptr[0]
        data = X["data"][a:ptr[-1]].astype(np.float32)
        ind = X["indices"][a:ptr[-1]].astype(np.int32)
        starts, ends = ptr[sel - sel[0]] - a, ptr[sel - sel[0] + 1] - a
        take = np.concatenate([np.arange(s, e) for s, e in zip(starts, ends)])
        datas.append(data[take])                       # copies, so the block can be freed
        inds.append(ind[take])
        lens.extend(ends - starts)
        del data, ind
        print(f"  rows {b0:,}-{min(b0 + block, n_obs):,}: {len(sel):,} selected", flush=True)
    indptr_new = np.concatenate([[0], np.cumsum(lens)])
    return sparse.csr_matrix((np.concatenate(datas), np.concatenate(inds), indptr_new), shape=(len(rows), n_var))


def gene_panel(var_names: pd.Index, hvg: np.ndarray, detect: np.ndarray, gmt: dict) -> list[str]:
    markers = yaml.safe_load((REPO / "metadata/annotation/markers.yaml").read_text())
    keep = set()
    for genes in markers["level1"].values():
        keep |= set(genes)
    for lin in markers["level2"].values():
        for genes in lin.values():
            keep |= set(genes)
    for genes in gmt.values():
        keep |= set(genes)
    em = REPO / "docs/phase1/sources/expression_map.tsv"
    if em.is_file():
        keep |= set(pd.read_csv(em, sep="\t")["gene"].astype(str))
    pos = pd.Series(np.arange(len(var_names)), index=var_names)
    core = [g for g in var_names if g in keep]
    hv = [var_names[i] for i in np.flatnonzero(hvg) if var_names[i] not in keep]
    hv = sorted(hv, key=lambda g: -detect[pos[g]])[:N_HVG]
    return [g for g in core + hv if detect[pos[g]] >= MIN_DETECT]


def encode_gene(col: np.ndarray, vmax: float) -> np.ndarray:
    """Bitmap (ceil(n/8) bytes, LSB first) + uint8 values of the nonzero cells, in cell order."""
    nz = col > 0
    bits = np.packbits(nz, bitorder="little")
    vals = np.clip(np.round(col[nz] / vmax * 255), 1, 255).astype(np.uint8)
    return np.concatenate([bits, vals])


def module_score(Xd: sparse.csc_matrix, genes: list[int], bins: np.ndarray, rng: np.random.Generator,
                 rows: np.ndarray) -> np.ndarray:
    """AddModuleScore-style score for the given rows: mean(set) - mean(expression-matched controls)."""
    ctrl = []
    for g in genes:
        pool = np.flatnonzero(bins == bins[g])
        ctrl.append(rng.choice(pool, min(50, len(pool)), replace=False))
    ctrl = np.unique(np.concatenate(ctrl))
    sub = Xd[rows]
    return np.asarray(sub[:, genes].mean(axis=1)).ravel() - np.asarray(sub[:, ctrl].mean(axis=1)).ravel()


# ── atlas build ──────────────────────────────────────────────────────────────
def donor_key_obs(f: h5py.File, rows: np.ndarray, rules: dict, skip=()) -> dict[str, pd.Series]:
    """obs columns that identify a cell's donor unit (pseudobulk_aggregate.cell_unit_keys), for the given rows."""
    cols = ("sample_id", "gse") + tuple(r["donor"]["column"] for r in rules["studies"].values()
                                        if r["donor"]["source"] == "obs")
    return {c: read_obs(f, c)[rows].reset_index(drop=True) for c in dict.fromkeys(cols)
            if c in f["obs"] and c not in skip}


def donor_table(atlas: str, obs: dict, rules: dict) -> pd.DataFrame:
    """One row per exported cell: its unit in the curated donor map (metadata/donor_map.csv)."""
    dm = pd.read_csv(REPO / "metadata/donor_map.csv", dtype=str)
    dm = dm[dm["atlas"] == atlas]
    dm["key"] = map_unit_keys(dm, rules).values
    keys = cell_unit_keys(pd.DataFrame(obs), rules)
    return dm.drop_duplicates("key").set_index("key").reindex(keys.values)


def cat_field(out: Path, fid: str, label: str, values, colours=None, group=None) -> dict:
    """Write obs/<fid>.bin (category codes, most frequent first) and return its meta.json field entry."""
    s = pd.Series(values).astype(str).replace({"nan": "NA", "None": "NA", "": "NA"})
    cats = s.value_counts().index.tolist()
    codes = pd.Categorical(s, categories=cats).codes
    dt = np.uint8 if len(cats) < 255 else np.uint16
    write(out / f"obs/{fid}.bin", codes.astype(dt))
    return {"id": fid, "label": label, "kind": "cat", "dtype": np.dtype(dt).name,
            "categories": cats, "counts": [int(x) for x in s.value_counts().values],
            "colours": colours(cats) if colours else palette(cats), "group": group}


def sample_rows(f: h5py.File, lab: pd.DataFrame, n_cells: int, rng: np.random.Generator) -> np.ndarray:
    """The exported subsample (level2 x study strata); the first draws of a generator seeded with SEED."""
    gse = read_obs(f, "gse").astype(str)
    strata = lab["level2"].astype(str).values + "|" + gse.values
    return subsample(pd.Series(strata), n_cells, rng)


def refresh_donor_fields(atlas: str) -> None:
    """Rewrite the donor-level obs fields (DONOR_FIELDS) of an existing build from the current donor map.

    The subsample is recomputed (same seed and strata) and checked against the exported Level 2 codes, so
    cell i stays the same cell; genes, scores and embeddings are untouched."""
    out = OUT / atlas.lower()
    meta = json.loads((out / "meta.json").read_text())
    rules = yaml.safe_load((REPO / "metadata/donor_rules.yaml").read_text())
    lab = pd.read_parquet(REPO / f"metadata/annotation/{atlas.lower()}_labels.parquet")
    with h5py.File(atlas_h5ad(atlas), "r") as f:
        rows = None
        for n_try in sorted({meta["n_cells"], 60_000, 150_000}):      # n_cells is the realised subsample size
            r = sample_rows(f, lab, n_try, np.random.default_rng(SEED))
            if len(r) == meta["n_cells"]:
                rows = r
                break
        if rows is None:
            raise SystemExit(f"{atlas}: could not reproduce the {meta['n_cells']:,}-cell subsample")
        obs = donor_key_obs(f, rows, rules)
    f2 = next(x for x in meta["fields"] if x["id"] == "level2")
    codes = np.frombuffer((out / "obs/level2.bin").read_bytes(), np.dtype(f2["dtype"]).newbyteorder("<"))
    if not np.array_equal(np.array(f2["categories"])[codes], lab["level2"].astype(str).values[rows]):
        raise SystemExit(f"{atlas}: recomputed subsample does not match the exported Level 2 labels")
    d = donor_table(atlas, obs, rules)
    pos = {x["id"]: i for i, x in enumerate(meta["fields"])}
    changed = []
    for fid, label in DONOR_FIELDS:
        new = cat_field(out, fid, label, d[fid].fillna("NA").values, None, "Donors (curated)")
        old = meta["fields"][pos[fid]]
        if (old["categories"], old["counts"]) != (new["categories"], new["counts"]):
            changed.append(fid)
        meta["fields"][pos[fid]] = new
    (out / "meta.json").write_text(json.dumps(meta, separators=(",", ":")))
    print(f"{atlas}: donor fields refreshed ({len(rows):,} cells); changed: {', '.join(changed) or 'none'}")


def build_atlas(atlas: str, n_cells: int) -> None:
    rng = np.random.default_rng(SEED)
    out = OUT / atlas.lower()
    if out.exists():
        shutil.rmtree(out)
    rules = yaml.safe_load((REPO / "metadata/donor_rules.yaml").read_text())
    lab = pd.read_parquet(REPO / f"metadata/annotation/{atlas.lower()}_labels.parquet")
    fields = COMMON + EXTRA[atlas]
    with h5py.File(atlas_h5ad(atlas), "r") as f:
        barcodes = read_obs(f, f["obs"].attrs["_index"]).astype(str).values
        gse = read_obs(f, "gse").astype(str)
        rows = sample_rows(f, lab, n_cells, rng)     # rng continues into the per-cell scores
        print(f"{atlas}: {len(rows):,} of {len(barcodes):,} cells", flush=True)
        obs = {c: read_obs(f, c)[rows].reset_index(drop=True) for _, _, c, _ in fields
               if c in f["obs"]}
        obs.update(donor_key_obs(f, rows, rules, skip=obs))
        obs["gse"] = gse[rows].reset_index(drop=True)
        umap = ad.io.read_elem(f["obsm"]["X_umap"])[rows]
        var_names = pd.Index(np.asarray(ad.io.read_elem(f["var"][f["var"].attrs["_index"]])).astype(str))
        hvg = np.asarray(ad.io.read_elem(f["var"]["highly_variable"])).astype(bool)
        X = read_rows(f, rows)
    n = len(rows)
    meta = {"atlas": atlas, "n_cells": n, "n_cells_total": int(len(barcodes)), "seed": SEED, "fields": []}

    # embeddings
    def put_xy(xy):
        lo, hi = np.nanmin(xy, 0), np.nanmax(xy, 0)
        q = np.round((xy - lo) / (hi - lo) * 65535).astype(np.uint16)
        return q, [float(lo[0]), float(hi[0]), float(lo[1]), float(hi[1])]
    q, rng_xy = put_xy(umap)
    write(out / "umap.bin", q.ravel())
    meta["umap_range"] = rng_xy
    meta["l2umap"] = {}
    bc = pd.Series(np.arange(n), index=barcodes[rows])
    for p in sorted((REPO / "results/annotation/level2").glob(f"{atlas.lower()}_*.parquet")):
        lin = p.stem[len(atlas) + 1:]
        u = pd.read_parquet(p, columns=["umap1", "umap2"])
        u = u[u.index.isin(bc.index)]
        ids = bc.loc[u.index].values.astype(np.uint32)
        order = np.argsort(ids)
        q, r = put_xy(u[["umap1", "umap2"]].values[order])
        write(out / f"l2umap/{lin}.bin", np.concatenate([ids[order].view(np.uint16), q.ravel()]))
        meta["l2umap"][lin] = {"n": int(len(ids)), "range": r}

    # labels
    def add_cat(fid, label, values, colours=None, group=None):
        meta["fields"].append(cat_field(out, fid, label, values, colours, group))

    def add_num(fid, label, values, group=None):
        v = pd.to_numeric(pd.Series(values), errors="coerce").values.astype(float)
        if not np.isfinite(v).any():
            return
        lo, hi = np.nanpercentile(v, [0.5, 99.5])
        write(out / f"obs/{fid}.bin", quantise(v, float(lo), float(hi)))
        meta["fields"].append({"id": fid, "label": label, "kind": "num", "dtype": "uint16",
                               "range": [float(lo), float(hi)], "group": group})

    l1 = lab["level1"].astype(str).values[rows]
    l2 = lab["level2"].astype(str).values[rows]
    add_cat("level1", "Level 1 (curated lineage)", l1,
            lambda cats: [LINEAGE_HUE.get(c, "#888888") for c in cats], "Cell types")
    l2_parent = pd.Series(l1).groupby(pd.Series(l2)).agg(lambda s: s.value_counts().index[0]).to_dict()

    def l2_colours(cats):
        col = {}
        for lin in set(l2_parent.values()):
            kids = sorted([c for c in cats if l2_parent.get(c) == lin])
            for c, h in zip(kids, hex_shades(LINEAGE_HUE.get(lin, "#888888"), len(kids))):
                col[c] = h
        return [col.get(c, "#888888") for c in cats]
    add_cat("level2", "Level 2 (curated sub-state)", l2, l2_colours, "Cell types")
    meta["level2_parent"] = l2_parent
    if atlas == "BM25":
        cnv = pd.read_parquet(REPO / "metadata/annotation/bm25_cnv_status.parquet")
        add_cat("cnv_status", "inferCNV call (caution: not specific)", cnv["cnv_status"].astype(str).values[rows],
                None, "Other")
        add_num("cnv_score", "inferCNV score", cnv["cnv_score"].values[rows], "Other")
    for fid, label, col, kind in fields:
        if col not in obs:
            continue
        grp = ("CellTypist" if fid.startswith("ct_") else "Quality control" if kind == "num" and fid != "age"
               else "Author labels" if fid in ("author_state", "author_annot", "author_type", "sctype")
               else "Samples and disease")
        if kind == "cat":
            add_cat(fid, label, obs[col].values, None, grp)
        else:
            add_num(fid, label, obs[col].values, grp)

    # donor-level metadata from the curated donor map
    d = donor_table(atlas, obs, rules)
    for fid, label in DONOR_FIELDS:
        add_cat(fid, label, d[fid].fillna("NA").values, None, "Donors (curated)")

    # genes
    Xc = X.tocsc()
    detect = np.diff(Xc.indptr) / n
    panel = gene_panel(var_names, hvg, detect, read_gmt(GMT))
    pos = pd.Series(np.arange(len(var_names)), index=var_names)
    genes_meta, total = [], 0
    budget = BUDGET_MB[atlas] * 1e6
    for g in panel:
        col = Xc[:, pos[g]].toarray().ravel()
        vmax = float(col.max())
        if vmax <= 0:
            continue
        sz = write(out / f"genes/{g}.bin.gz", encode_gene(col, vmax), gz=True)
        total += sz
        genes_meta.append([g, round(vmax, 4), round(float(detect[pos[g]]), 4)])
        if total > budget:
            print(f"  gene budget reached at {len(genes_meta)} genes", flush=True)
            break
    meta["genes"] = genes_meta
    meta["gene_bytes"] = total
    print(f"  {len(genes_meta)} genes, {total / 1e6:.1f} MB", flush=True)

    # per-cell cytokine-response scores
    mapping = yaml.safe_load((REPO / "metadata/celltype_mapping/hucira_level2.yaml").read_text())["celltype_mapping"]
    gmt = read_gmt(GMT)
    mods = read_modules()
    mod_rank = mods.drop_duplicates("module").set_index("module")["module_rank"]
    mean = np.asarray(X.mean(axis=0)).ravel()
    expressed = mean > 0
    ranks = pd.Series(mean[expressed]).rank(method="first").values
    bins = np.full(len(mean), -1)
    bins[expressed] = np.minimum((ranks - 1) * 25 // expressed.sum(), 24).astype(int)
    dct = pd.Series(l2).map(mapping)
    cytokines, modules = {}, {}
    for ct in sorted(dct.dropna().unique()):
        cells = np.flatnonzero((dct == ct).values)
        if len(cells) < 20:
            continue
        for name, genes in gmt.items():
            c_t, cyto = name.split("|", 1)
            if c_t != ct:
                continue
            gi = [pos[g] for g in genes if g in pos.index and expressed[pos[g]]]
            if len(gi) < 5:
                continue
            s = module_score(Xc, gi, bins, rng, cells)
            cytokines.setdefault(cyto, {})[ct] = (cells, s)
        for mod, g in mods[mods["celltype"] == ct].groupby("module"):
            genes = set().union(*[set(gmt.get(f"{ct}|{c}", [])) for c in g["cytokine"]])
            gi = [pos[x] for x in genes if x in pos.index and expressed[pos[x]]]
            if len(gi) >= 5:
                modules.setdefault(mod, {})[ct] = (cells, module_score(Xc, gi, bins, rng, cells))

    def put_scores(kind, table, order=None):
        lst = []
        for name, per_ct in sorted(table.items(), key=lambda kv: (order(kv[0]) if order else 0, kv[0])):
            z = np.full(n, -127, np.int8)
            for ct, (cells, s) in per_ct.items():
                sd = s.std() or 1.0
                z[cells] = np.clip(np.round((s - s.mean()) / sd * 25), -126, 127).astype(np.int8)
            fn = f"{kind}__{name}".replace("/", "_").replace(" ", "_")
            write(out / f"scores/{fn}.bin.gz", z, gz=True)
            lst.append({"name": name, "file": fn, "celltypes": sorted(per_ct)})
        return lst
    meta["cytokine_scores"] = put_scores("cytokine", cytokines)
    # modules in display order (consensus groups); "grouped" = one of the huCIRA highly similar groups
    meta["module_scores"] = [dict(r, grouped=r["name"] in GROUP_LABEL.values()) for r in
                             put_scores("module", modules, lambda m: int(mod_rank.get(m, 999)))]
    meta["score_scale"] = 25        # int8 / 25 = z within cell type; -127 = not scored
    meta["dictionary_celltypes"] = sorted({ct for v in cytokines.values() for ct in v})
    meta["level2_to_dictionary"] = {k: v for k, v in mapping.items() if k in set(l2)}
    meta["markers"] = yaml.safe_load((REPO / "metadata/annotation/markers.yaml").read_text())
    (out / "meta.json").write_text(json.dumps(meta, separators=(",", ":")))
    size = sum(p.stat().st_size for p in out.rglob("*") if p.is_file())
    print(f"{atlas}: wrote {out} ({size / 1e6:.1f} MB)")


# ── shared results ───────────────────────────────────────────────────────────
def build_results() -> None:
    out = OUT / "results"
    out.mkdir(parents=True, exist_ok=True)

    def dump(name, df, cols, rnd=None):
        df = df[[c for c in cols if c in df.columns]].copy()
        for c in df.select_dtypes("float").columns:
            df[c] = df[c].round((rnd or {}).get(c, 4))
        recs = json.loads(df.to_json(orient="split", index=False))
        (out / f"{name}.json").write_text(json.dumps({"columns": recs["columns"], "data": recs["data"]},
                                                    separators=(",", ":")))
    rd = lambda p: pd.read_csv(REPO / p, sep="\t")  # noqa: E731
    # tested: the cell-count rule (scripts/test_eligibility.py); untested rows keep their effect, q = NaN
    dump("meta", tested(rd("docs/phase1/meta_main.tsv")),
         ["group", "celltype", "cytokine", "k", "studies", "n_test", "n_ref", "effect", "ci_low", "ci_high", "p", "q",
          "i2", "tested"])
    dump("study_effects", tested(rd("docs/phase1/study_effects_main.tsv")),
         ["group", "celltype", "cytokine", "study", "condition", "timepoint", "n_test", "n_ref", "median_cells_test",
          "median_cells_ref", "effect", "se", "score_p", "tested"], rnd={"median_cells_test": 1, "median_cells_ref": 1})
    mods = read_modules()
    ranks = mods.drop_duplicates("module")[["module", "module_rank", "module_grouped"]]
    dump("module_meta", tested(rd("docs/phase2/modules/module_meta.tsv")).merge(ranks, on="module", how="left"),
         ["group", "celltype", "module", "module_rank", "module_grouped", "k", "studies", "n_test", "n_ref", "effect",
          "ci_low", "ci_high", "p", "q", "tested"])
    dump("propeller", rd("docs/phase2/composition/propeller.tsv"),
         ["contrast", "compartment", "study", "condition", "celltype", "prop_test", "prop_ref", "log2_ratio", "p",
          "q", "n_test", "n_ref"])
    studies = studies_current_donors(rd("docs/data_collection/atlas_studies.tsv"))
    dump("studies", studies, list(studies.columns))
    dump("modules", mods, ["celltype", "cytokine", "module", "module_rank", "module_grouped", "n_genes",
                           "frac_ifn_shared", "frac_cell_cycle"])
    build_ch_cell(out)
    print(f"wrote {out}")


def tested(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise the `tested` flag to bool and enforce q = NaN for untested rows (they never enter BH)."""
    if "tested" not in df.columns:
        raise SystemExit("result table has no `tested` column: re-run the donor-level analysis")
    df = df.copy()
    df["tested"] = df["tested"].map({True: True, False: False, "True": True, "False": False}).astype(bool)
    if "q" in df.columns:
        df.loc[~df["tested"], "q"] = np.nan
    return df


def studies_current_donors(st: pd.DataFrame) -> pd.DataFrame:
    """atlas_studies.tsv with donors_by_condition recomputed from the current metadata/donor_map.csv
    (as scripts/summarise_data_collection.py does), so a relabelled or excluded donor shows at once;
    a study's entry is replaced only when its counts differ (keeps the original order otherwise)."""
    dm = pd.read_csv(REPO / "metadata/donor_map.csv", dtype=str)
    dm["condition"] = dm["condition"].fillna("unassigned")
    dm = dm[dm["donor"].notna()].drop_duplicates(["atlas", "study", "donor"])
    parse = lambda s: dict((k, int(v)) for k, v in (x.rsplit(" ", 1) for x in str(s).split(", ") if " " in x))  # noqa: E731
    st = st.copy()
    for i, r in st.iterrows():
        vc = dm.loc[(dm["atlas"] == r["atlas"]) & (dm["study"] == r["study"]), "condition"].value_counts()
        if len(vc) and parse(r["donors_by_condition"]) != vc.to_dict():
            st.at[i, "donors_by_condition"] = ", ".join(f"{k} {v}" for k, v in vc.items())
    return st


def build_ch_cell(out: Path) -> None:
    """Original cell-level huCIRA results for the Heimlich 2024 CH cohort, and their donor-level counterparts.

    Cytokines are labelled with the project's response modules (docs/phase1/response_modules.tsv): the huCIRA
    paper's highly similar cytokine groups (Supplementary Table 11), fixed across cell types, other cytokines
    singly. A signature not in that table (e.g. too few genes for the donor-level tests) gets its module by the
    same fixed rule. module_rank is the display order (response_modules.module_order).
    """
    mods = read_modules()[["celltype", "cytokine", "module"]]
    runs = {"all": REPO / "results/hucira_pbmc4_gse/CXG_HEIMLICH_2024",
            "women": REPO / "results/hucira_pbmc4_gse_by_sex/female/CXG_HEIMLICH_2024",
            "men": REPO / "results/hucira_pbmc4_gse_by_sex/male/CXG_HEIMLICH_2024"}
    rows = []
    for sex, root in runs.items():
        for f in sorted(root.glob("*_vs_control/*/tables/*_all_enrichment_results.csv")):
            geno = f.parts[-4].split("_vs_")[0]
            d = pd.read_csv(f)
            rows.append(pd.DataFrame({
                "contrast": f"CH_{geno} vs healthy", "sex": sex, "lfc": d["threshold_lfc"], "direction": d["direction"],
                "celltype": d["celltype_signature"], "cytokine": d["cytokine"],
                "NES": d["NES"], "p": d["NOM p-val"],
                "q": d["FDR q-val"], "n_cells_test": d["num_cells_1"], "n_cells_ref": d["num_cells_2"],
                "n_genes": d["num_shared_genes_signature"],
                "lead": d["Lead_genes"].fillna("").map(lambda s: ";".join(str(s).split(";")[:25]))}))
    cell = pd.concat(rows, ignore_index=True)
    # the original tables hold repeated runs of some tests; keep the first, as in hucira_concordance.read_original
    cell = cell.drop_duplicates(["contrast", "sex", "lfc", "celltype", "cytokine"])
    cell = cell.merge(mods, on=["celltype", "cytokine"], how="left")
    cell["module"] = cell["module"].fillna(module_of(cell["cytokine"]))
    rank = {m: i for i, m in enumerate(module_order(cell[["cytokine", "module"]].drop_duplicates()))}
    cell["module_rank"] = cell["module"].map(rank)
    cell["module_grouped"] = cell["module"].isin(set(GROUP_LABEL.values()))
    dump_table(out / "ch_cell.json", cell, rnd=6)
    conc = []
    for name, lfc in (("main", 1.0), ("lfc0.8", 0.8)):
        meta = tested(pd.read_csv(REPO / f"docs/phase1/meta_{name}.tsv", sep="\t"))[
            ["group", "celltype", "cytokine", "q", "ci_low", "ci_high", "tested"]]
        c = pd.read_csv(REPO / f"docs/phase2/hucira_prerank/concordance_{name}.tsv", sep="\t")
        c = c[c["contrast"].astype(str).str.startswith("CXG_HEIMLICH_2024__")].copy()
        c["group"] = c["contrast"].str.split("__").str[1].str.replace("_vs_", " vs ")
        c = c.merge(meta.rename(columns={"q": "q_donor_score", "tested": "tested_donor"}),
                    on=["group", "celltype", "cytokine"], how="left")
        c["tested_donor"] = c["tested_donor"].eq(True)   # no donor-level row: not tested
        conc.append(c.assign(lfc=lfc)[["group", "lfc", "celltype", "cytokine", "NES_cell", "q_cell", "n_test", "n_ref",
                                       "g_donor", "ci_low", "ci_high", "score_p", "q_donor_score",
                                       "tested_donor"]].rename(columns={"group": "contrast"}))
    dump_table(out / "ch_concordance.json", pd.concat(conc, ignore_index=True), rnd=6)
    print(f"CH cell level: {len(cell)} tests; concordance rows: {sum(len(c) for c in conc)}")


def dump_table(path: Path, df: pd.DataFrame, rnd: int = 4) -> None:
    df = df.copy()
    for c in df.select_dtypes("float").columns:
        df[c] = df[c].round(rnd)
    recs = json.loads(df.to_json(orient="split", index=False))
    path.write_text(json.dumps({"columns": recs["columns"], "data": recs["data"]}, separators=(",", ":")))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--atlas", choices=["BM25", "PBMC4"])
    ap.add_argument("--n-cells", type=int, default=150_000)
    ap.add_argument("--results-only", action="store_true")
    ap.add_argument("--donor-fields", choices=["BM25", "PBMC4"],
                    help="only refresh the donor-level fields of an existing atlas build from metadata/donor_map.csv")
    a = ap.parse_args()
    if a.donor_fields:
        refresh_donor_fields(a.donor_fields)
        return
    if a.results_only or not a.atlas:
        build_results()
    if a.atlas:
        build_atlas(a.atlas, a.n_cells)


if __name__ == "__main__":
    main()
