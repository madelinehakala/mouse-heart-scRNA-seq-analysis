library(dplyr)
library(Seurat)
library(patchwork)

zone2.data = Read10X(data.dir = "/Users/madelinehakala/mouse-heart-scRNA-seq-analysis/dataZone2/")
zone3.data = Read10X(data.dir = "/Users/madelinehakala/mouse-heart-scRNA-seq-analysis/dataZone3/")

zone2 = CreateSeuratObject(counts = zone2.data, project = "Zone II", min.cells = 3, min.features = 200)
zone2

zone3 = CreateSeuratObject(counts = zone3.data, project = "Zone III", min.cells = 3, min.features = 200)
zone3
