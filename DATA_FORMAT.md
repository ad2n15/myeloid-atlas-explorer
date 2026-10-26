# Explorer data format

Written by `scripts/build_explorer.py`. All binaries are little-endian. `n` = number of
subsampled cells (`meta.json: n_cells`); cell `i` is the same cell in every file.

| File | Content |
|---|---|
| `data/<atlas>/meta.json` | `fields` (id, label, kind `cat`/`num`, dtype, categories, counts, colours, range, group), `genes` ([name, max, fraction detected]), `cytokine_scores` / `module_scores` ([name, file, celltypes]), `umap_range`, `l2umap`, `level2_parent`, `level2_to_dictionary`, `markers`, `score_scale` |
| `umap.bin` | uint16 x,y interleaved (2n values); value / 65535 maps onto `umap_range` = [xmin, xmax, ymin, ymax] |
| `l2umap/<lineage>.bin` | uint32 cell ids (m values), then uint16 x,y interleaved for those cells (range in `l2umap[lineage].range`) |
| `obs/<field>.bin` | `cat`: category codes (uint8, or uint16 if ≥ 255 categories). `num`: uint16, 0 = missing, 1..65535 = `range[0]..range[1]` (0.5–99.5th percentiles; values outside are clipped) |
| `genes/<gene>.bin.gz` | gzip of: nonzero bitmap (ceil(n/8) bytes, bit i = cell i, least-significant bit first), then one uint8 per nonzero cell in cell order; value = byte / 255 × gene max (log1p of counts per 10,000) |
| `scores/<kind>__<name>.bin.gz` | gzip of int8 per cell: −127 = not scored (cell type has no such signature); otherwise byte / `score_scale` = z-score within the cell type |
| `data/results/*.json` | donor-level result tables as `{columns, data}` (row arrays) |
| `data/results/ch_cell.json` | original cell-level huCIRA results (gseapy prerank) for the Heimlich 2024 CH cohort: contrast, sex (all/women/men), signature threshold (lfc), cell type, cytokine, response module of that cytokine in that cell type (`docs/phase1/response_modules.tsv`; "unassigned" if the signature has no module there), NES, p, FDR q, cells compared, leading genes (first 25); duplicated runs removed |
| `data/results/ch_concordance.json` | the same tests matched to the donor-level analysis: cell NES and q next to donor Hedges' g with 95% CI, score-test p and meta-analysis q |

Per-cell scores: mean log1p of the signature genes minus the mean of expression-matched
control genes (25 bins of mean expression, 50 controls per signature gene), computed for cells
whose curated Level 2 maps to the huCIRA dictionary cell type
(`metadata/celltype_mapping/hucira_level2.yaml`), then z-scored within that cell type.
