library(dplyr)
library(Seurat)
library(patchwork)

## Allowing RScript to access variables from Python script
args = commandArgs(trailingOnly = TRUE)
directory = args[1]
sampleLabel = args[2]
markerGene1 = args[3]
markerGene2 = args[4]

## Path to data
data = Read10X(data.dir = directory)

## Analysis
sample = CreateSeuratObject(counts = data, project = sampleLabel, min.cells = 3, min.features = 200)
sample
sample[["percent.mt"]] = PercentageFeatureSet(sample, pattern = "^MT-")
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 = FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
sample = subset(sample, subset = nFeature_RNA > 1000 & nFeature_RNA < 4500)
sample =  NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)
if (is.null(markerGene2)) {
  VlnPlot(sample, features = markerGene1)
} else {
  VlnPlot(sample, features = c(markerGene1, markerGene2))
}
sampleExpressionData = GetAssayData(object = sample)
sampleExpressionDf = as.data.frame(sampleExpressionData)
sampleExpressionDf = t(sampleExpressionDf)
rownames(sampleExpressionDf) = colnames(sample)
write.csv(sampleExpressionDf, file = sprintf('%s_expression_data.csv', sampleLabel))
