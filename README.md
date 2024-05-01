# Building a Pipeline for scRNA-Seq Analysis

## Overview
The overall goal for this project was to build a pipeline that assesses the level of gene expression heterogeneity exhibited by the muscle tissue of the mouse heart (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132658). In striving to complete this goal, however, we built a readily generalizable pipeline. As long as data is in the correct format, this pipeline can be used to analyze a wide variety of scRNA-Seq datasets. 
**For detailed project design information, see the Design Document in the Wiki tab.**

## The Pipeline
1. Retrieve data.
2. Assess expression levels of known marker genes for your cell type of interest (for us, these were Myh6 and Actc1, which are marker genes for muscle cells).
3. Select for your cell type of interest, filtering out other cell types, by utilizing the expression levels of the marker genes.
4. Calculate the average expression levels of genes that you know to be implicated in your cell/tissue type of interest.
5. Calculate the relative expression levels of the above genes in individual cells when compared to the population average.
6. Output these findings to a well-organized csv file for additional analysis.
7. Data visualization and additional analysis.

## Required Dependencies
- Seurat (R)
- patchwork (R)
- Pandas (Python)
- Seaborn (Python)

## Running Our Code
To run our code, utilize the following command:
```
python3 runExperiment.py -d [DATADIRECTORY] -g [GENESOFINTERESTFILE] -s [SAMPLELABEL] -l [LOGFILENAME] -g1 [MARKERGENE1] -g2 [MARKERGENE2] -f [FILTEREDOUTFILENAME] -p [PLOTFILENAME] -n [NORMALIZEDOUTFILENAME]
```
- [DATADIRECTORY]: directory where your data is stored
- [GENESOFINTERESTFILE]: text file that contains your genes of interest that you would like to include in the final output
- [SAMPLELABEL]: label that describes your data (ex: Zone2MouseHeart)
- [LOGFILENAME]: the name you would like to give your logfile
- [MARKERGENE1]: a gene that is known to be a "marker" for your cell type of interest
- [MARKERGENE2]: a second gene that is known to be a "marker" for your cell type of interest (this argument is optional)
- [FILTEREDOUTFILENAME]: the name you would like to give the file that contains absolute gene expression data for your cell type of interest
- [PLOTFILENAME]: the name you would like to give the outputted plots
- [NORMALIZEDOUTFILENAME]: the name you would like to give the file that contains normalized gene expression data for your cell type of interest

Ex: When we ran this script using data from Zone II of the mouse heart, we used the following command (all data and scripts were stored in the indicated directory):
```
python3 runExperiment.py -d dataZone2 -g muscle_contraction_genes.txt -s Zone2 -l Zone2LogFile.txt -g1 Myh6 -g2 Actc1 -f zone2_muscle_cells.csv -p zone2_muscle_cell_expression_data.png -n normalized_zone2_muscle_cells.csv
```
## Data Formatting
For Seurat (and the rest of our code) to work, data must be in a very specific format and placed in the corresponding directories. 

Ex: One of data directories:
dataZone2 -> barcodes.tsv.gz features.tsv.gz matrix.mtx.gz

To download the data that we analyzed for this project in the required format, run the following command:
```
sh getData.sh
```
If you are using your own data, just make sure it follows the appropriate format.

## Tutorial
Follow the tutorial below (using the tutorial data found in this repository) to practice working with this pipeline.

Step 1:
```
git clone https://github.com/madelinehakala/mouse-heart-scRNA-seq-analysis
```

Step 2:
```
cd mouse-heart-scRNA-seq-analysis
```

Step 3:
```
python3 runExperiment.py -d tutorial_data -g tutorial_genes_of_interest.txt -s tutorial -l tutorialLogFile.txt -g1 Myh6 -g2 Actc1 -f tutorial_muscle_cells.csv -p tutorial_muscle_cell_expression_data.png -n tutorial_normalized_muscle_cells.csv
```
**When prompted to enter expression cutoffs for Myh6 and Actc1, enter 0.75 and 3, respectively.**
