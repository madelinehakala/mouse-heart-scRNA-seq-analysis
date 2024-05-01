# Building a Pipeline for scRNA-Seq Analysis

## Overview
This project utilizes scRNA-Seq data from the mouse heart, which was collected by the Barefield Lab (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132658). Our overall goal is to build a pipeline that assesses the level of gene expression homogeneity exhibited by the muscle tissue of the heart. **For detailed project design information, see the Design Document in the Wiki tab.**

## The Pipeline
1. Retrieve data from GEO.
2. Assess expression levels of known muscle marker genes (Myh6 and Actc1).
3. Select for muscle cells, filtering out other cell types, by utilizing the expression levels of the muscle marker genes.
4. Calculate the average expression levels of 20-30 muscle contraction genes (specifically in the muscle cells).
5. Calculate the relative expression levels of the 20-30 muscle contraction genes in individual cells compared to the population average.
6. Output these findings to a well-organized csv file for additional analysis.
7. Data visualization and additional analysis.

## Data Download
For Seurat (and the rest of our code) to work, data must be in a very specific format and placed in the corresponding directories like so:

dataZone2 -> barcodes.tsv.gz features.tsv.gz matrix.mtx.gz

dataZone3 -> barcodes.tsv.gz features.tsv.gz matrix.mtx.gz

To download the data yourself in the required format, run the following command:
```
sh getData.sh
```

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
