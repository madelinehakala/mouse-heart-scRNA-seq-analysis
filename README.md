# Assessing Gene Expression Homogeneity in Mouse Heart Sarcomeres: scRNA-Seq Analysis

# Overview
This project utilizes scRNA-Seq data from the mouse heart, which was collected by the Barefield Lab (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132658). Our overall goal is to build a pipeline that assesses the level of homogeneity exhibited by the muscle tissue of the heart. **For detailed project design information, see the Design Document in the Wiki tab.**

# The Pipeline
1. Retrieve data from GEO.
2. Quantify gene expression levels.
3. Select for muscle cells, filtering out other cell types, by utilizing the expression levels of muscle marker genes (Myh6 and Actc1).
4. Calculate the average expression levels of 20-30 muscle contraction genes (specifically in the muscle cells).
5. Calculate the relative expression levels of the 20-30 muscle contraction genes in individual cells compared to the population average.
6. Output these findings to a well-organized table for additional analysis.

# Running Our Code
To run our code, utilize the following command:
```
python3 wrapper.py -o [OUTPUT DIRECTORY] -l [LOG FILE NAME]
```
# Required Dependencies
LIST REQUIRED PACKAGES HERE
