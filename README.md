# Demonstration: Pathway redistribution reveals a shared signaling backbone and context-dependent regulatory modules in RNA-binding protein networks

---

This repository demonstrates a deep learning approach for predicting gene expression levels using nucleic acid-binding site data (derived from ChIP-seq) and gene co-expression networks, as described in our study:

**“Pathway redistribution reveals a shared signaling backbone and context-dependent regulatory modules in RNA-binding protein networks.”**
[bioRxiv](https://doi.org/10.1101/2025.03.03.641203)

In this analysis, DNA- or RNA-binding sites derived from co-expression data are used as inputs to the deep learning model.
- DNA-binding site data are used to predict gene expression levels in four cell types including **HepG2**.
- RNA-binding site data are used to predict gene expression levels in four cell types including **K562**.

This work builds upon our previous study:

**"Systematic discovery of directional regulatory motifs associated with human insulator sites."**
[bioRxiv](https://doi.org/10.1101/2024.01.20.573595)
[GitHub](https://github.com/reposit2/insulator)

Due to file size limitations on GitHub, the required dataset is hosted on Zenodo.
[Zenodo repository](https://zenodo.org/records/15281537).

---

## Deep Learning Analysis — Predicting Gene Expression Levels

### Requirements

#### Data Files
Download the dataset from [Zenodo](https://zenodo.org/records/15661042):

- `data.tgz`

Then extract it using:

```bash
tar zxf data.tgz
```
     
Create output directories: 
     
```bash
mkdir -p train_out/3celltypes_hepg2
mkdir -p train_out/3celltypes_k562      
```

### Execution

Run the main script:

```bash
python deeplearning.py  
```     

### Outputs

- Spearman’s correlation coefficients between actual and predicted gene expression levels
- DeepLIFT contribution scores for nucleic acid-binding sites

All output files are saved in the `./train_out/` directory.  
A complete archive of the results (`train_out.tgz`) is available from [Zenodo](https://zenodo.org/records/15661042).

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

## Identification of Putative Regulatory Target Genes

Create an output directory:

```bash
mkdir geneidlist
```

To generate gene-level rankings based on DeepLIFT scores, edit the following scripts as instructed in the embedded comments and run:

```bash
perl dlscore_ranksort_with_ensembl_csv.pl
perl rankscore.pl
```

The DeepLIFT score files used in this analysis are included in the `./train_out/3celltypes_k562_original` directory within the `train_out.tgz` archive.

### Outputs

Ensembl gene IDs representing putative regulatory targets for each NABP, identified based on normalized DeepLIFT score rankings.  
The results are saved in the `./geneidlist/` directory. A subset of these results is also included in the `geneidlist` directory of this repository.

---

## Functional Analysis: Ontology & Semantic Inference with ChatGPT

### A. Functional Enrichment Analysis (PANTHER Overrepresentation Test)

[PANTHER](https://pantherdb.org/)

#### Steps:

1. Upload your gene list (e.g., `K562_PKM_rank2_2024-10-08_03-21-22_332_gid.txt`).
2. Configure the following options:
   - **Organism:** `Homo sapiens`
   - **Analysis Type:** `Statistical overrepresentation test`
   - **Annotation Data Sets:**
     - PANTHER Protein Class
     - PANTHER GO-Slim Biological Process
     - PANTHER GO-Slim Cellular Component
     - PANTHER GO-Slim Molecular Function
     - Reactome pathways
3. Upload a **reference list**, such as `fourcelltypeswithhffk562_fpkm_gid.txt` or `fourcelltypeswithhffhepg2_fpkm_gid.txt`, available in this repository.
4. Select:
   - **Test Type:** Fisher’s Exact
   - **Correction:** False Discovery Rate (FDR)
5. Click **Launch analysis**
6. Export the result as a `Table`.

---

### B. Semantic Function Prediction Using ChatGPT-4o

#### Step 1: Convert Ensembl IDs to Gene Symbols

Edit the file paths for `$infile` and `$outfile` in the script below, then run:

```bash
mkdir chatgpt
cd chatgpt
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
perl genelist.pl
```

#### Step 2: Predict Functions via ChatGPT

Use the gene list in the prompt template:
`protein_system_analysis_instructions_gptinput2_nonprotein.txt`

Reference: Engineered prompt from Hu et al., *Nature Methods* 2024
[Link to paper](https://doi.org/10.1038/s41592-024-02525-x)

#### Step 3: Individual Gene Analysis

Use prompt template:
`gene_set_analysis_instructions.txt`
(Input 1–3 filled in as instructed)

#### Step 4: Literature Retrieval Using GPT-4 API

- Copy the ChatGPT output to a Word file to preserve formatting.
- Convert it to Markdown:

```bash
pandoc -s ChatGPT_genefunction.docx --wrap=none --extract-media=pandoc_out -t gfm -o ChatGPT_genefunction.md
```

#### Step 5: Generate Summary Table

Edit the file paths for `$infile` and `$outfile` in the script below, then run the script to parse the Markdown output and create a gene-function table:

```bash
perl gpt2table.pl
```

Then move the file:

```bash
mkdir data
cp -p ChatGPT_genefunction_table.tsv data/omics_revamped_LLM_DF_test.tsv
```

#### Environment Setup for GPT-4 API Access

Refer to the [Code Ocean README](https://codeocean.com/capsule/2519655/tree/v1) for setup instructions.

---

## LLM-Based Reference Search

Set up environment:

```bash
conda create -n llm_eval python=3.11.5
conda activate llm_eval
conda env config vars set OPENAI_API_KEY="<your api key>"
conda deactivate
conda activate llm_eval
echo $OPENAI_API_KEY  # confirm key setup
```

```python
# Python test
import os
import openai
openai.api_key = os.environ["OPENAI_API_KEY"]
```

Clone and install:

```bash
git clone git@github.com:idekerlab/llm_evaluation_for_gene_set_interpretation.git
cd llm_evaluation_for_gene_set_interpretation
pip install -r requirements.txt
```


Edit the configuration file `jsonFiles/reference_checking_revision_test.json` and update the file paths for `LLM_analysisFilePath` and `toSaveFilePath` in the script below.  
Then, run the script to perform the reference search:

```bash
python 4reference_search_and_validation_test.py
less data/omics_revamped_LLM_ref_DF_test.tsv
```

---

## Functional Analysis: Gene Set Enrichment Analysis (GSEA) Using DeepLIFT Scores

We provide an example of Gene Set Enrichment Analysis (GSEA) using DeepLIFT-derived gene rankings to quantify pathway-level redistribution of regulatory influence.

This example corresponds to the analysis shown in **Figure 6**, **Table 4** and **Supplementary Table S7-S17**.

### Overview

Genes are ranked based on PKM-associated DeepLIFT scores comparing **K562** and **NPC** cells.  
GSEA is then applied to identify pathways enriched toward either side of the ranked distribution.

Unlike overrepresentation analysis, this approach uses the full ranked gene list and does not require an arbitrary cutoff.

### Input

- DeepLIFT score matrices (gene-wise, normalized):
  - `DNA_3_norm_ensembl.csv` (K562)
  - `DNA_2_norm_ensembl.csv` (NPC)
- Orientation map:
  - `orient_map.csv`
- Background gene list:
  - `fourcelltypeswithhffk562_fpkm_id.csv`
- Gene set library (Reactome):
  - `reactome/reactome_v2026_ensembl.gmt`

### Example command

Set up environment:
```bash
pip install gseapy
```
Run GSEA:
```bash
python gsea_rownorm_one_factor_bg_posneg.py \
  --factor PKM \
  --matrix_A ./3celltypes_k562/2024-10-08_03-21-22/DeepLIFT/DNA_3_norm_ensembl.csv \
  --matrix_B ./3celltypes_k562/2024-10-08_03-21-22/DeepLIFT/DNA_2_norm_ensembl.csv \
  --label_A K562 \
  --label_B NPC \
  --orient_map orient_map.csv \
  --background fourcelltypeswithhffk562_fpkm_id.csv \
  --library reactome/reactome_v2026_ensembl.gmt \
  --min_size 10 \
  --max_size 5000 \
  --min_nonzero 3 \
  --nperm 10000 \
  --threads 7 \
  --weight 1 \
  --rownorm minmax \
  --ranking_mode raw \
  --jitter 1e-12 \
  --save_leading_edge \
  --outdir out_PKM_rowminmax_GSEA_reactomev2026_ensembl
```

---

## Reference and License Information

### DEcode  
Tasaki, S., et al.
**Deep learning decodes the principles of differential gene expression.**
*Nature Machine Intelligence* (2020).
[Link](https://doi.org/10.1038/s42256-020-0201-6)[PubMed](https://pubmed.ncbi.nlm.nih.gov/32671330/) — BSD 3-Clause License

### DeepLIFT  
Shrikumar, A., et al.
**Learning Important Features Through Propagating Activation Differences.**
*arXiv* (2017).
[Link](https://arxiv.org/abs/1704.02685) — MIT License

### PANTHER  
Thomas, P.D., et al.
**PANTHER: Making genome-scale phylogenetics accessible to all.**
*Protein Science* (2022).
[Link](https://doi.org/10.1002/pro.4218) — [Website](https://pantherdb.org/)

### Language Models  
Hu, M., et al.
**Evaluation of large language models for discovery of gene set function.**
*Nature Methods* (2024).
[Link](https://doi.org/10.1038/s41592-024-02525-x)
[GitHub](https://github.com/idekerlab/llm_evaluation_for_gene_set_interpretation)
[Code Ocean](https://doi.org/10.24433/CO.7045777.v1) — MIT License
