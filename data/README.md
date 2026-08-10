# Data files

All four `reference_sets_*.csv` files and `chilas_block_inventory.csv` are
**transcribed from the manuscript**, not measured here:

| file | manuscript source |
|---|---|
| `reference_sets_manual.csv`    | Table 4, "Manual Mapping" column |
| `reference_sets_qfacet_fm.csv` | Table 4, "Fast Marching (qFacet FM)" column |
| `reference_sets_qfacet_kd.csv` | Table 4, "K-Dimensional tree (qFacet KD tree)" column |
| `reference_sets_dse.csv`       | Table 4, "Discontinuity Set Extractor (DSE)" column |
| `chilas_block_inventory.csv`   | Table 5, block ID / source sector / volume |

They are provided so that the validation and block-statistics stages are
runnable and auditable. They are **not** a substitute for the primary data:

* the 65-million-point Chilas virtual outcrop model,
* the 4,818 manually mapped discontinuity orientations (Section 4.2),
* the 1,472 DSE and 1,872 / 16,500 qFacet plane orientations (Section 4.3),
* the CloudCompare 2.5D volume projects behind Table 5, and
* the RocFall3 model and its outputs (Section 4.5).

None of those were supplied with the code archive. Until they are, Sections
4.1 to 4.3 and 4.5 of the manuscript remain unverifiable from the repository.
