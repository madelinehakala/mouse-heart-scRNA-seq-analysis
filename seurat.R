library(dplyr)
library(Seurat)
library(patchwork)

## Paths to data
zone2.data = Read10X(data.dir = "dataZone2/")
zone3.data = Read10X(data.dir = "dataZone3/")

## Zone 2 Analysis
zone2 = CreateSeuratObject(counts = zone2.data, project = "Zone II", min.cells = 3, min.features = 200)
zone2
zone2[["percent.mt"]] = PercentageFeatureSet(zone2, pattern = "^MT-")
VlnPlot(zone2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 = FeatureScatter(zone2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(zone2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
zone2 = subset(zone2, subset = nFeature_RNA > 1000 & nFeature_RNA < 4500)
zone2 =  NormalizeData(zone2, normalization.method = "LogNormalize", scale.factor = 10000)
VlnPlot(zone2, features = c("Myh6", "Actc1"))
zone2ExpressionData = GetAssayData(object = zone2)
zone2ExpressionDf = as.data.frame(zone2ExpressionData)
zone2ExpressionDf = t(zone2ExpressionDf)
rownames(zone2ExpressionDf) = colnames(zone2)
write.csv(zone2ExpressionDf, file = 'zone2_expression_data.csv')

## Zone 3 Analysis
zone3 = CreateSeuratObject(counts = zone3.data, project = "Zone III", min.cells = 3, min.features = 200)
zone3
zone3[["percent.mt"]] = PercentageFeatureSet(zone3, pattern = "^MT-")
VlnPlot(zone3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 = FeatureScatter(zone3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(zone3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
zone3 = subset(zone3, subset = nFeature_RNA > 1000 & nFeature_RNA < 4500)
zone3 =  NormalizeData(zone3, normalization.method = "LogNormalize", scale.factor = 10000)
VlnPlot(zone3, features = c("Myh6", "Actc1"))
zone3ExpressionData = GetAssayData(object = zone3)
zone3ExpressionDf = as.data.frame(zone3ExpressionData)
zone3ExpressionDf = t(zone3ExpressionDf)
rownames(zone3ExpressionDf) = colnames(zone3)
write.csv(zone3ExpressionDf, file = 'zone3_expression_data.csv')


