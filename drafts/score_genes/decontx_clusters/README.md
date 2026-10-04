# DecontX cluster labels

The population-group labels supplied to DecontX as its `z` argument, one file per
slice. They are an input to the ambient-RNA correction, not an output, and are
not regenerated deterministically by a re-run, so they are kept here to make the
correction reproducible exactly.

| File | Slice | Cells | Clusters |
|---|---|---|---|
| `slice_1_clusters.csv` | 1 (L321, tumor-bearing) | 124,956 | 25 |
| `slice_2_clusters.csv` | 2 (L321, tumor-bearing) | 126,485 | 25 |
| `slice_3_clusters.csv` | 3 (L321, control) | 63,997 | 25 |
| `slice_4_clusters.csv` | 4 (L34, control) | 57,086 | 25 |
| `slice_5_clusters.csv` | 5 (L34, tumor-bearing) | 206,210 | 25 |
| `slice_6_clusters.csv` | 6 (L34, tumor-bearing) | 267,373 | 25 |

All six slices use the same partition: library-size normalization and log1p,
30 truncated-SVD components, then MiniBatchKMeans with k = 25 (seed 0). A Leiden
partition was not used because building the neighbor graph exceeds available
memory on the larger slices, and one method is applied to every slice.

Each file has a single column, `cluster`. Row *i* is the label of the *i*-th cell
of that slice's exported count matrix, which contains the cells retained after
low-transcript QC in their original order; the files carry no cell identifiers.

`sensitivity/slice_1_leiden5_clusters.csv` is the earlier five-cluster Leiden
partition of slice 1 (Scanpy; 30 PCs, 15 neighbors, resolution 0.5). It is not used
for the correction itself, only as one of the alternatives compared by
`thesis_plots/decontx_partition_sensitivity.py`.

The partitions are written by `score_genes/write_decontx_clusters.py`.
`score_genes/run_decontx_correct.py` exports the count matrices and
`score_genes/run_decontx.R` reads these labels when running DecontX.
