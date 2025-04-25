# Demonstration: Inferring the cell type–specific functional roles of nucleic acid-binding proteins using deep learning on co-expression networks

---

This repository provides a demonstration of predicting gene expression levels using DNA/RNA-binding sites derived from ChIP-seq data and gene co-expression data, as described in our study:

**“Inferring the cell type–specific functional roles of nucleic acid-binding proteins using deep learning on co-expression networks.”**
[bioRxiv](https://doi.org/10.1101/2025.03.03.641203)

In this analysis, the input to the deep learning model includes DNA- or RNA-binding sites derived from gene co-expression data.  
- DNA-binding site data are used to predict gene expression levels in four cell types including **HepG2**.  
- RNA-binding site data are used to predict gene expression levels in four cell types including **K562**.

This deep learning approach builds upon our previous study:

**"Systematic discovery of directional regulatory motifs associated with human insulator sites."**
[bioRxiv](https://doi.org/10.1101/2024.01.20.573595)
[GitHub](https://github.com/reposit2/insulator)

Due to file size limitations on GitHub, the required dataset is hosted on Zenodo.
[Zenodo repository](https://zenodo.org/records/15281537).

---

### Requirements

#### Data Files
Please download the following file from the [Zenodo repository](https://doi.org/10.5281/zenodo.8216164)

- `data.tgz`

Extract it with the following command:

```bash
tar zxf data.tgz
```

### Run

To execute the analysis, run the following script:

```bash
python deeplearning.py
```

### Outputs
- Spearman’s correlation coefficients between actual and predicted gene expression levels
- DeepLIFT contribution scores for DNA/RNA-binding sites involved in the prediction

These outputs are saved in the `./train_out/` directory.

```
cat ./train_out/3celltypes_hepg2/2025-04-25_14-43-19/test_data/cor_tbl.txt 
sample  estimate        statistic       p.value method  alternative
HepG2   0.798046041634558       502914646.020844        0       Spearman's rank correlation rho two.sided
HFF     0.797082158458548       505314950.378294        0       Spearman's rank correlation rho two.sided
HMEC    0.80094973464168        495683741.746181        0       Spearman's rank correlation rho two.sided
MedianExp       0.842166455760396       393044046.642756        0       Spearman's rank correlation rho two.sided
NPC     0.802257203936158       492427824.096745        0       Spearman's rank correlation rho two.sided

cat ./train_out/3celltypes_k562/2025-04-25_18-00-57/test_data/cor_tbl.txt
sample  estimate        statistic       p.value method  alternative
HFF     0.808626750957343       476566097.436869        0       Spearman's rank correlation rho two.sided
HMEC    0.808591626333202       476653566.323638        0       Spearman's rank correlation rho two.sided
K562    0.830602496151738       421841128.414551        0       Spearman's rank correlation rho two.sided
MedianExp       0.852274296673504       367873055.809033        0       Spearman's rank correlation rho two.sided
NPC     0.811765877663049       468748905.791841        0       Spearman's rank correlation rho two.sided
```

---

### Reference and License

#### DEcode
Tasaki, S., Gaiteri, C., Mostafavi, S., and Wang, Y.  
**Deep learning decodes the principles of differential gene expression.**  
*Nature Machine Intelligence* (2020).  
[Link to paper](https://doi.org/10.1038/s42256-020-0201-6)  
License: BSD 3-Clause

#### DeepLIFT
Shrikumar, A., Greenside, P., and Kundaje, A.  
**Learning Important Features Through Propagating Activation Differences.**  
*arXiv* (2017).  
[Link to paper](https://arxiv.org/abs/1704.02685)
License: MIT

