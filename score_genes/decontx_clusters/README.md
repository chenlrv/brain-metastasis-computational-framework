# DecontX cluster labels

The population-group labels supplied to DecontX as its `z` argument, one file per
slice. They are an input to the ambient-RNA correction, not an output, and are
not regenerated deterministically by a re-run, so they are kept here to make the
correction reproducible exactly.

| File | Slice | Cells | Clusters | Method |
|---|---|---|---|---|
| `slice_1_clusters.csv` | 1 (L321, tumor-bearing) | 124,956 | 5 | Leiden (Scanpy; 30 PCs, 15 neighbors, resolution 0.5) |
| `slice_2_clusters.csv` | 2 (L321, tumor-bearing) | 126,485 | 25 | k-means on 30 truncated-SVD components |
| `slice_3_clusters.csv` | 3 (L321, control) | 63,997 | 6 | Leiden (as slice 1) |
| `slice_4_clusters.csv` | 4 (L34, control) | 57,086 | 25 | k-means (as slice 2) |
| `slice_5_clusters.csv` | 5 (L34, tumor-bearing) | 206,210 | 25 | k-means (as slice 2) |
| `slice_6_clusters.csv` | 6 (L34, tumor-bearing) | 267,373 | 25 | k-means (as slice 2) |

Each file has a single column, `cluster`. Row *i* is the label of the *i*-th cell
of that slice's exported count matrix, which contains the cells retained after
low-transcript QC in their original order; the files carry no cell identifiers.

The k-means partitions are written by `score_genes/write_decontx_clusters.py`.
`score_genes/run_decontx_correct.py` exports the count matrices and
`score_genes/run_decontx.R` reads these labels when running DecontX.
