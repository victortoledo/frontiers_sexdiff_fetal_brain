# Codes and files for the original research article "Sex Differences in Gene Regulatory Networks During Mid-Gestational Brain Development"

## Data folders:

### raw data:
contains the 120 mid-gestational brain homogenate samples analyzed throughout this article
data was obtained in the repository https://figshare.com/articles/Summary_statistics_for_expression_quantitative_trait_loci_in_the_developing_human_brain_and_their_enrichment_in_neuropsychiatric_disorders/6881825

### metadata:
contains the covariates associated to the samples, obtained during e-mail exchange with corresponding author

## Analysis:

### 0. Preprocessing:
contains the preparation of the expression data input to build co-expression and regulatory networks
- checks for data stratification by covariate (sex, age, RIN and batch), filters low variable genes, discards outliers, corrects of batch effect

### 1. CEMiTool:
creates and analyzes CEMiTool co-expression modules
- builds the co-expression network, calculates class specific activity, runs enrichment of biological processes in each of the modules, integrates PPI data

### 2. netZoo: 
builds and analyzes gene regulatory networks

#### 2.1. PANDA
builds the networks integrating RNA counts, PPI and _motif_ data

#### 2.2. CONDOR + ALPACA
creates differential modules and extracts biological meaning of differences in interactions
- CONDOR: clusters PANDA network elements into modules
- ALPACA: compares community structure between two modules
- enrichment of modules in biological processes databases, extraction of transcription factors and genes that drive these differences, enrichment of modules into cell type, disorders and hormone targets databases
