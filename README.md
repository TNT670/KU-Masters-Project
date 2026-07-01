# Cancer Tissue Microarray Classification

This project applies machine learning classifiers (decision tree, SVM, MLP, etc.) to Affymetrix microarray gene expression data to distinguish cancerous tissue samples from healthy tissue samples. Three tissue types are covered: breast, lung, and colon. Data was downloaded from the [NIH Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/) database.

---

## Data sources

Each GEO Series (GSE) entry represents one batch of samples, which may contain healthy tissue, cancerous tissue, or a mix. Data files are distributed as tarballs of Affymetrix `.CEL` files. All datasets use the **Affymetrix HG-U133 Plus 2.0** array platform (GEO accession: GPL570).

| Cancer type | GSE datasets |
|---|---|
| Breast | [GSE5460](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5460), [GSE10780](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10780), [GSE26457](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26457), [GSE29431](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29431), [GSE42568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568), [GSE66162](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66162) |
| Lung   | [GSE19188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19188), [GSE19804](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804), [GSE27262](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27262), [GSE51024](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51024) |
| Colon  | [GSE23878](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23878), [GSE71571](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71571), [GSE92921](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92921), [GSE143985](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143985) |

---

## Directory structure

```
./
├── preprocess.R            # Parameterized R preprocessing pipeline (run this)
├── config.yaml             # Dataset definitions, sample exclusions, output names
├── classifier_utils.py     # Shared classifier helpers (config/label loading, feature pruning, plots)
├── counts.txt              # Reference: raw sample counts per tarball (pre-exclusion)
├── exclusions.txt          # Reference: per-GSE sample exclusion lists
├── breast_cancer/
│   ├── GSE*_RAW.tar         # one raw tarball per GSE
│   ├── breast_cancer.ipynb  # classifier notebook (complete + leave-one-out modes)
│   ├── breast_labels.json   # per-GSE label config
│   └── *_txt.txt            # text-label files for text-parsed GSEs
├── lung_cancer/
│   ├── GSE*_RAW.tar
│   ├── lung_cancer.ipynb
│   ├── lung_labels.json
│   └── *_txt.txt
└── colon_cancer/
    ├── GSE*_RAW.tar
    ├── colon_cancer.ipynb
    └── colon_labels.json
```

Tar files are extracted automatically into a `.extracted/` subdirectory within each cancer directory on the first run. Delete the tar files or directories to force re-extraction.

---

## Dependencies

### R

Install Bioconductor packages (required for preprocessing):

```r
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("affy", "frma", "biomaRt", "WGCNA", "preprocessCore", "sva", "hgu133plus2frmavecs"))
```

SVA may require you to install `genefilter`, `limma`, and `edgeR`:

```r
BiocManager::install(c("genefilter", "limma", "edgeR"))
```

`optparse` and `yaml` (required by `preprocess.R`) are installed automatically if missing.

For offline probe-to-gene mapping via the `--local-annot` flag (optional):
```r
BiocManager::install("hgu133plus2.db")
```

Also install for exploratory analysis:
```r
install.packages(c("ggplot2", "ggfortify"))
```

### Python

A virtual environment is recommended. conda can help simplify GPU setup for TensorFlow. Otherwise pip, uv, etc. are acceptable alternatives.

Required packages:
- `numpy`, `pandas`, `matplotlib`
- `scikit-learn`
- `keras` (Keras 3) + `keras-tuner`, and at least one Keras backend (`torch`, `tensorflow`, or `jax`)
- `xgboost`

Optional but useful:
- `sklearnex` — accelerates certain scikit-learn operations on Intel hardware
- `joblib` — cache trained models to skip retraining on re-runs

#### Keras backend

Keras 3 can use a backend of your choice among PyTorch, TensorFlow, or JAX, and you must have one installed. Each classifier notebook selects the backend near the top of the imports cell:

```python
os.environ['KERAS_BACKEND'] = 'torch'
```

Change `torch` to `tensorflow` or `jax` to match your install, or remove the line and instead export `KERAS_BACKEND` in your shell before launching Jupyter. The notebooks default to `torch`.

If you are using legacy Keras 2 (where Keras is just `tf.keras` on TensorFlow), this line has no effect and can be left in or removed. Just ensure TensorFlow is installed in that case. The model-building code targets the Keras 3 API.

---

## How to run

**Summary of the run order:**

1. Generate the `.dat` table(s) with `preprocess.R` for the mode you want.
2. Open the corresponding Jupyter notebook.
3. Set `HELD_OUT_GSE`.
4. Run all cells.

### Step 1 — R preprocessing

**All-inclusive run** (combine all GSEs into one training matrix):

```bash
Rscript preprocess.R --cancer breast_cancer
Rscript preprocess.R --cancer lung_cancer
Rscript preprocess.R --cancer colon_cancer
```

Output: `{cancer_dir}/{output_prefix}.dat`

**Leave-one-out (LOO) run** (hold out one GSE as a test set):

```bash
# Hold out GSE42568 from breast cancer
Rscript preprocess.R --cancer breast_cancer --loo 42568
```

Output:
- `{output_prefix}_loo{gse}_train.dat` — training samples (all GSEs except held-out)
- `{output_prefix}_loo{gse}_test.dat` — test samples (held-out GSE only, batch-corrected to match training scale)

**Local annotation (offline, does not query Ensembl):**

```bash
Rscript preprocess.R --cancer breast_cancer --local-annot
```

Uses the bundled `hgu133plus2.db` package instead of querying Ensembl. This option is included in the case that Ensembl is unavailable to connect. Although faster, mappings may differ slightly from the biomaRt output due to different curation pipelines.

**Custom config path:**

```bash
Rscript preprocess.R --cancer breast_cancer --config /path/to/config.yaml
```

### Step 2 — Python classifiers

The classifier phase reads the `.dat` tables from Step 1, selects a subset informative genes, and trains and evaluates seven classifiers (decision tree, random forest, SVM, XGBoost, KNN, MLP, logistic regression). There is one notebook per cancer type, each notebook using the shared logic file `classifier_utils.py`. Each notebook documents every stage inline.

```
breast_cancer/breast_cancer.ipynb   +  breast_cancer/breast_labels.json
lung_cancer/lung_cancer.ipynb       +  lung_cancer/lung_labels.json
colon_cancer/colon_cancer.ipynb     +  colon_cancer/colon_labels.json
classifier_utils.py                 (shared: config loading, label building, feature pruning, plots)
```

**Choose the run mode** with the single `HELD_OUT_GSE` variable in the first code cell:

- `HELD_OUT_GSE = None` for a **complete run**. Reads `{prefix}.dat`, combines all GSEs into a single dataset, and reports 5-fold cross-validated accuracy per model. Run `preprocess.R --cancer <type>` first.
- `HELD_OUT_GSE = '42568'` (a GSE number as a string) for a **leave-one-out run.** Reads `{prefix}_loo{GSE}_train.dat` and `{prefix}_loo{GSE}_test.dat`, trains on the other GSEs, and tests on the held-out GSE. Run `preprocess.R --cancer <type> --loo <GSE>` first to generate the matching tables. Leave-one-out runs cache fitted models under `pkls/` so re-runs are fast.

> The MLP uses Keras 3. Make sure a backend is installed and selected (see [Keras backend](#keras-backend) above).

---

## Output format

Each `.dat` file is a space-delimited table with:
- Rows: individual tissue samples
- Columns: Entrez gene IDs (one column per gene) + a `Batch` column indexed by GSE
- Values: log-scale, batch-corrected, quantile-normalized expression intensities

---

## What the preprocessing does

Raw Affymetrix `.CEL` files contain fluorescence intensities for tens of thousands of probe sequences on a chip. Each probe is designed to hybridize with a specific RNA transcript, so its intensity reflects how actively that gene is being expressed in the sample. The pipeline converts these raw intensities into a clean, normalized, gene-level expression matrix through five steps:

1. **Frozen RMA normalization (fRMA).** Standard RMA normalization adjusts intensities relative to all other arrays in the same batch, which means arrays from different studies can't be compared directly since each was normalized against a different set. Frozen RMA (fRMA) solves this by normalizing every array against a single pre-computed reference distribution derived from a large external collection. This enables combining samples across many studies.

2. **Probe ID to Entrez gene ID mapping.** Probe IDs are Affymetrix-specific and not directly useful for biology. biomaRt queries the Ensembl database to map each probe ID to its corresponding Entrez gene ID. Probes with no known gene mapping are discarded.

3. **Collapsing duplicate probes per gene.** Multiple probes on the array may map to the same gene. WGCNA's `collapseRows` selects a single representative probe per gene (the one with the highest mean expression across all samples) and discards the rest, producing one value per gene per sample.

4. **Quantile normalization.** This forces all samples to share the same empirical expression distribution, removing residual intensity-level differences between arrays that aren't attributable to biology.

5. **ComBat batch effect removal.** Even after the above steps, samples from different labs can still have technical differences in reagent lots, scanner calibration, or sample handling, collectively referred to as batch effects. ComBat models these effects and removes them while preserving biological variation. Each GSE is treated as one batch.

For LOO runs, quantile normalization is performed separately between the LOO GSE and the collective remaining GSEs. In addition, ComBat batch effect removal is performed first on the training GSEs, then separately on the test GSE; data from the test GSE is only referenced against the training GSEs when removing batch effects from the test GSE. This prevents any information from the test set from influencing the correction.

---

## Adding a new dataset

1. Download the `_RAW.tar` file from GEO and place it in the appropriate cancer directory. The series **must** use the Affymetrix HG-U133 Plus 2.0 array platform (GEO accession GPL570), the same platform as the other datasets. A platform mismatch will likely cause errors in the preprocessing script, since the fRMA reference and probe-to-gene mappings are platform-specific.
2. Add an entry to `config.yaml` under the relevant cancer type with the GSE number and any samples to exclude.
3. Re-run `preprocess.R` for that cancer type.

---

## References

This project builds on the following open-source libraries:

### R / Bioconductor (preprocessing)

- **affy** — Gautier L, Cope L, Bolstad BM, Irizarry RA (2004). affy—analysis of Affymetrix GeneChip data at the probe level. *Bioinformatics* 20(3):307–315. [doi:10.1093/bioinformatics/btg405](https://doi.org/10.1093/bioinformatics/btg405)
- **frma (frozen RMA)** — McCall MN, Bolstad BM, Irizarry RA (2010). Frozen robust multiarray analysis (fRMA). *Biostatistics* 11(2):242–253. [doi:10.1093/biostatistics/kxp059](https://doi.org/10.1093/biostatistics/kxp059)
  - frmaTools / `hgu133plus2frmavecs` — McCall MN, Irizarry RA (2011). Thawing Frozen Robust Multi-array Analysis (fRMA). *BMC Bioinformatics* 12:369. [doi:10.1186/1471-2105-12-369](https://doi.org/10.1186/1471-2105-12-369)
- **biomaRt** — Durinck S, Spellman PT, Birney E, Huber W (2009). Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. *Nature Protocols* 4(8):1184–1191. [doi:10.1038/nprot.2009.97](https://doi.org/10.1038/nprot.2009.97)
- **WGCNA** — Langfelder P, Horvath S (2008). WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics* 9:559. [doi:10.1186/1471-2105-9-559](https://doi.org/10.1186/1471-2105-9-559)
  - `collapseRows` — Miller JA, Cai C, Langfelder P, et al. (2011). Strategies for aggregating gene expression data: the collapseRows R function. *BMC Bioinformatics* 12:322. [doi:10.1186/1471-2105-12-322](https://doi.org/10.1186/1471-2105-12-322)
- **preprocessCore (quantile normalization)** — Bolstad BM, Irizarry RA, Åstrand M, Speed TP (2003). A comparison of normalization methods for high density oligonucleotide array data based on variance and bias. *Bioinformatics* 19(2):185–193. [doi:10.1093/bioinformatics/19.2.185](https://doi.org/10.1093/bioinformatics/19.2.185)
- **sva (ComBat)** — Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD (2012). The sva package for removing batch effects and other unwanted variation in high-throughput experiments. *Bioinformatics* 28(6):882–883. [doi:10.1093/bioinformatics/bts034](https://doi.org/10.1093/bioinformatics/bts034)
  - ComBat method — Johnson WE, Li C, Rabinovic A (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics* 8(1):118–127. [doi:10.1093/biostatistics/kxj037](https://doi.org/10.1093/biostatistics/kxj037)
- **hgu133plus2.db** — Carlson M. *hgu133plus2.db: Affymetrix Human Genome U133 Plus 2.0 Array annotation data.* Bioconductor annotation package.
- **ggplot2** — Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis.* Springer-Verlag New York. <https://ggplot2.tidyverse.org>
- **ggfortify** — Tang Y, Horikoshi M, Li W (2016). ggfortify: Unified Interface to Visualize Statistical Results of Popular R Packages. *The R Journal* 8(2):474–485. [doi:10.32614/RJ-2016-060](https://doi.org/10.32614/RJ-2016-060)

### Python (classifiers)

- **scikit-learn** — Pedregosa F, et al. (2011). Scikit-learn: Machine Learning in Python. *Journal of Machine Learning Research* 12:2825–2830. [link](https://jmlr.org/papers/v12/pedregosa11a.html)
- **NumPy** — Harris CR, et al. (2020). Array programming with NumPy. *Nature* 585:357–362. [doi:10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)
- **pandas** — McKinney W (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 56–61. [doi:10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a)
- **Matplotlib** — Hunter JD (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering* 9(3):90–95. [doi:10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55)
- **XGBoost** — Chen T, Guestrin C (2016). XGBoost: A Scalable Tree Boosting System. *Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining*, 785–794. [doi:10.1145/2939672.2939785](https://doi.org/10.1145/2939672.2939785)
- **Keras** — Chollet F, et al. (2015). Keras. <https://keras.io>
- **KerasTuner** — O'Malley T, Bursztein E, Long J, Chollet F, Jin H, Invernizzi L, et al. (2019). KerasTuner. <https://github.com/keras-team/keras-tuner>
- **PyTorch** (Keras 3 backend) — Paszke A, et al. (2019). PyTorch: An Imperative Style, High-Performance Deep Learning Library. *Advances in Neural Information Processing Systems* 32:8024–8035. [link](https://papers.neurips.cc/paper/9015-pytorch-an-imperative-style-high-performance-deep-learning-library)
