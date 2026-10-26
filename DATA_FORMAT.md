# Explorer data format

Written by `scripts/build_explorer.py`. All binaries are little-endian. `n` = number of
subsampled cells (`meta.json: n_cells`); cell `i` is the same cell in every file.

| File | Content |
|---|---|
| `data/<atlas>/meta.json` | `fields` (id, label, kind `cat`/`num`, dtype, categories, counts, colours, range, group), `genes` ([name, max, fraction detected]), `cytokine_scores` ([name, file, celltypes]) / `module_scores` ([name, file, celltypes, grouped]; in display order, `grouped` = a huCIRA highly similar group), `umap_range`, `l2umap`, `level2_parent`, `level2_to_dictionary`, `markers`, `score_scale` |
| `umap.bin` | uint16 x,y interleaved (2n values); value / 65535 maps onto `umap_range` = [xmin, xmax, ymin, ymax] |
| `l2umap/<lineage>.bin` | uint32 cell ids (m values), then uint16 x,y interleaved for those cells (range in `l2umap[lineage].range`) |
| `obs/<field>.bin` | `cat`: category codes (uint8, or uint16 if ≥ 255 categories). `num`: uint16, 0 = missing, 1..65535 = `range[0]..range[1]` (0.5–99.5th percentiles; values outside are clipped) |
| `genes/<gene>.bin.gz` | gzip of: nonzero bitmap (ceil(n/8) bytes, bit i = cell i, least-significant bit first), then one uint8 per nonzero cell in cell order; value = byte / 255 × gene max (log1p of counts per 10,000) |
| `scores/<kind>__<name>.bin.gz` | gzip of int8 per cell: −127 = not scored (cell type has no such signature); otherwise byte / `score_scale` = z-score within the cell type |
| `data/results/*.json` | donor-level result tables as `{columns, data}` (row arrays); `module_meta` and `modules` carry `module_rank` (display order: the huCIRA consensus groups, `response_modules.module_order`) and `module_grouped`; `studies` is `docs/data_collection/atlas_studies.tsv` with `donors_by_condition` recomputed from `metadata/donor_map.csv` |
| `data/results/meta.json`, `module_meta.json`, `study_effects.json` | boolean `tested` (cell-count rule, `scripts/test_eligibility.py`): a contrast × cell type is reported and enters the BH correction only if the donor pseudobulks hold a median of ≥ `min_median_cells_tested` (`pseudobulk_config.yaml`, 100) cells in each arm (test and reference). Untested rows keep their effect and CI, with `q` = null. `study_effects` also carries `median_cells_test` / `median_cells_ref` (median cells per donor pseudobulk in each arm); a tested meta-analysis pools only its tested studies |
| `data/results/ch_cell.json` | original cell-level huCIRA results (gseapy prerank) for the Heimlich 2024 CH cohort: contrast, sex (all/women/men), signature threshold (lfc), cell type, cytokine, response module of that cytokine (`docs/phase1/response_modules.tsv`; huCIRA highly similar groups, fixed across cell types, other cytokines singly; signatures not in that table get their module by the same rule), `module_rank`, `module_grouped`, NES, p, FDR q, cells compared, leading genes (first 25); duplicated runs removed |
| `data/results/ch_concordance.json` | the same tests matched to the donor-level analysis: cell NES and q next to donor Hedges' g with 95% CI, score-test p, meta-analysis q (`q_donor_score`) and `tested_donor` (the donor-level cell-count rule; false also when there is no donor-level row; `q_donor_score` is null when false) |

Per-cell scores: mean log1p of the signature genes minus the mean of expression-matched
control genes (25 bins of mean expression, 50 controls per signature gene), computed for cells
whose curated Level 2 maps to the huCIRA dictionary cell type
(`metadata/celltype_mapping/hucira_level2.yaml`), then z-scored within that cell type.

Response modules: the huCIRA paper's highly similar cytokine groups (Oesinghaus et al. 2025,
Supplementary Table 11; `metadata/hucira_signatures/hucira_cytokine_groups.tsv`), the same in
every cell type; a module's per-cell score uses the union of its members' signature genes in
that cell type. Every other cytokine is its own module, named by its display name.

Donor fields (`donor`, `condition`, `stage`, … in `obs/`) come from `metadata/donor_map.csv`; after a
change there, `python3 scripts/build_explorer.py --donor-fields <ATLAS>` rewrites them for an existing
build (the subsample is recomputed with the same seed and checked against the exported Level 2 labels).
The GSE285943 donor labelled 'cytopenia (MF?)' has condition `cytopenia_MF_uncertain` and is excluded
from the donor-level analyses (`metadata/donor_rules.yaml: exclude_conditions`); its cells stay in the maps.
