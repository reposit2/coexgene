# Demonstration: Prediction of gene expression levels based on DNA/RNA-binding sites

---

This repository provides a demonstration of predicting gene expression levels using DNA/RNA-binding sites derived from ChIP-seq data and gene co-expression data, as described in our study:

**“Functional prediction of DNA/RNA-binding proteins using deep learning based on gene expression correlation.”**
[bioRxiv](https://doi.org/10.1101/2025.03.03.641203)

In this analysis, the input to the deep learning model includes DNA- or RNA-binding sites derived from gene co-expression data.  
- DNA-binding site data are used to predict gene expression levels in four cell types including **HepG2**.  
- RNA-binding site data are used to predict gene expression levels in four cell types including **K562**.

This deep learning approach builds upon our previous study:

**"Systematic discovery of directional regulatory motifs associated with human insulator sites."**
[bioRxiv](https://doi.org/10.1101/2024.01.20.573595)
[GitHub](https://github.com/reposit2/insulator)

---

### Run

To execute the analysis, run the following script:

```bash
python deeplearning.py
```

### Outputs
- Spearman’s correlation coefficients between actual and predicted gene expression levels
- DeepLIFT contribution scores for DNA/RNA-binding sites involved in the prediction
These outputs are saved in the `./train_out/` directory.

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

