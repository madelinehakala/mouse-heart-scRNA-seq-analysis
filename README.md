# Building a Pipeline for scRNA-Seq Analysis

## Overview
Tje overall goal for this project is to build a pipeline that assesses the level of gene expression homogeneity exhibited by the muscle tissue of the heart (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132658). In striving to complete this goal, however, we built a pipeine that is generalizable beyond the muscle cells of the heart. As long as data is in the correct format, this pipeline can be used to analyze a wide variety of scRNA-Seq datasets. 
**For detailed project design information, see the Design Document in the Wiki tab.**

## The Pipeline
1. Retrieve data.
2. Assess expression levels of known marker genes for your cell type of interest (for us, these were Myh6 and Actc1).
3. Select for your cell type of interest, filtering out other cell types, by utilizing the expression levels of the marker genes.
4. Calculate the average expression levels of genes that you know to be implicated in your cell/tissue type of interest.
5. Calculate the relative expression levels of the above genes in individual cells compared to the population average.
6. Output these findings to a well-organized csv file for additional analysis.
7. Data visualization and additional analysis.

## Data Download
For Seurat (and the rest of our code) to work, data must be in a very specific format and placed in the corresponding directories. Example:

dataZone2 -> barcodes.tsv.gz features.tsv.gz matrix.mtx.gz

To download the data that we analyzed for this project in the required format, run the following command:
```
sh getData.sh
```
If you are using your own data, just make sure it follows the appropriate format.

## Running Our Code
To run our code, clone this repo, download data and required dependencies, and then utilize the following command:
```
python3 runExperiment.py -d [DIRECTORY] -l [LOG FILE NAME]
```
Ex: When I ran this script, I used the following command (all data and scripts were stored in the indicated directory):
```
python3 runExperiment.py -d /Users/madelinehakala/mouse-heart-scRNA-seq-analysis -l logFile.txt
```

## Required Dependencies
- Seurat (R)
- patchwork (R)
- Pandas (Python)
- Seaborn (Python)
