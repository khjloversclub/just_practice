### learn about the process how to make Rds files
#[BRIC Bio통신원] [당신의 논문 동료] scRNA-seq data 분석법 – Plot 4 종류부터 DEG까지,
# https://www.ibric.org/s.do?CopuYdakft

## 1. prologue (skip)
## 2. download scRNA-seq data & set up
# raw data are in the link above

## 2.4 make Seurat object
#  install tidyverse; skip
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

## 3. QC
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#nFeature_RNA이 200 초과, 2500 미만, 미토콘드리아 비율이 5% 미만인 세포만을 선별
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## Normalize & DimPlot
# make DimPlot
pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc)

pbmc <- ScaleData(pbmc)

pbmc <- RunPCA(pbmc) # << dimension is 50(seems the cause of error)

pbmc <- FindNeighbors(pbmc, dims = 1:10) 
# 본 튜토리얼에서는 PCA 1~10차원 정보를 사용함. 데이터에 따라 20차원 또는 50차원 정보까지 사용하는 경우가 있고, 주로 30차원 정도를 쓰는 것 같습니다. 튜토리얼 홈페이지에서 차원 결정 기준을 잘 참고해주시길 바랍니다.

pbmc <- FindClusters(pbmc, resolution = 0.5) 
# 본 튜토리얼에서는 resolution parameter로 0.5를 사용했습니다. Resolution parameter는 클러스터 수를 결정하는 인자임. 사용자가 임의로 결정해주는 것이며, 이 parameter를 결정하는 절대적인 기준은 없습니다. 논문에 따라 0.005부터 1.2까지 사용하는 것을 봤습니다. 

pbmc <- RunUMAP(pbmc, dims = 1:10) ###ERROR###
# 위의 FindNeighbors에서 사용한 차원가 동일한 차원을 사용하시면 됩니다.
# ------------------------------------------------
# x2set(Xsub, n_neighbors, metric, nn_method = nn_sub, n_trees, 에서 다음과 같은 에러가 발생했습니다: Non-finite entries in the input matrix
pbmc <- ScaleData(pbmc)  # 데이터 정규화
pbmc <- RunPCA(pbmc, npcs = 10)  # PCA 수행, origin: pbmc <- RunPCA(pbmc)/ '경고: Number of dimensions changing from 50 to 10' is printed out
pbmc <- FindNeighbors(pbmc, dims = 1:10)  # UMAP 이전에 필수
pbmc <- RunUMAP(pbmc, dims = 1:10)  # UMAP 실행
### error resolved
# ------------------------------------------------

DimPlot(pbmc, reduction = "umap")

## 4. save Seurat file (and be a Taljuninja)
setwd("~/Rstd_khj/GitHub_khj/just_practice/pbmc_folder") # save at the high rank file of WD(filtered_gene~~)
saveRDS(pbmc, "pbmc_tutorial.rds")

## load Seurat file(~.rds) & make 4 kinds of plots

# load Seurat file; just click 'pbmc_tutorial.rds' and OK

# 4.2. draw DimPlot, FeaturePlot, VlnPlot, DotPlot 
library(Seurat)

# 1) Dimplot
DimPlot(pbmc_tutorial)
DimPlot(pbmc_tutorial, pt.size = 1) # in/decrease the size of point 
DimPlot(pbmc_tutorial, cols = c("grey", "grey", "grey", "green", "grey", "grey", "grey", "grey", "grey")) #change the color !!! the colums start from 0(zero) !!!
DimPlot(pbmc_tutorial, label = T) # write the label the TOP of the UMAP
# Erase the UGLY x, y axes and legend
DimPlot(pbmc_tutorial) + NoAxes() + NoLegend()

# 2) draw FeatureFlot ; examine the amount of gene exp on UMAP
FeaturePlot(pbmc_tutorial,"SREBF1")
# order(bring) the buried cells gene SREBF1 expressed on the TOP and inc the size
FeaturePlot(pbmc_tutorial, "SREBF1", pt.size = 2, order = T) #...no B(bottom)
# then we can locate the gene expressing SRE~

# let's look for the other gene MS4A1, the B lympocyte marker
FeaturePlot(pbmc_tutorial, 'MS4A1')
# well... It's expressed significantly on the cluster 3(turn back the previous page if you forget the cluster numbers)

# investigate the amnt of exp of 2 gene at the same time
FeaturePlot(pbmc_tutorial, c('MS4A1','CD79A'),blend = T) # DO NOT FORGET to attach the term c, column!!!!!!!!!!

# 3) draw Violin plot = VlnPlot
VlnPlot(pbmc_tutorial, "MS4A1")
# wanna erase all points, make the size zero
VlnPlot(pbmc_tutorial, "MS4A1",pt.size = 0)

# EXTRA stage: draw BoxPlot
# the package Seurat do not offer BoxPlot option so we should load 'ggplot2'
library(ggplot2)
VlnPlot(pbmc_tutorial, "MS4A1",pt.size = 0) + geom_boxplot()

# 4) draw DotPlot; estimate the amount of various genes at the same time
# and SELECT what are the MARKERS available for each genes
DotPlot(pbmc_tutorial, features = c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")) + RotatedAxis()
# !!! DotPlot must contain 'RotatedAxis()'
# Change the color
DotPlot(pbmc_tutorial, features = c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4"), cols = c("lightgrey", "red")) + RotatedAxis()
# reverse the X, Y axis ; coord_flip()
DotPlot(pbmc_tutorial, features = c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4"), cols = c("lightgrey", "red")) + RotatedAxis() + coord_flip()

## 5. Cluster ANNOTATION
# First, identify the relation of each cluster >> drawing DimPlot again is the best&fast way.
DimPlot(pbmc_tutorial,label = T, pt.size = 0.5)
# find the similarity of clusters ; 0246 / 157 / 3 / 8                                                       (BTW, the cluster 0246 looks like a chicken drumstick)
# next, figure out the difference of 4 cluster
# i) extract the marker no.3 and no.8 clusters and do Annotation
# ii) compare no.0246 cluster with the rest clusters to find the features the 0246 cluster (and do the same at 157 clstr)
# iii) compare and find features of their subgroups within cluster 0246 and 157

# look over the marker genes already known via FeatureFlot
# it's easy to predict the rough functions of clusters before checking out some each type of marker of celltype
# we use the marker based of scRNA-seq offered the link below
# https://panglaodb.se/markers.html?cell_type=%27choose%27
# Now we browse some markers regarding Blood & Immune System (cause we're targeting PBMCs like T-B cells, NK cells! :D)
'''
T cell: THEMIS, CD3E, IL7R …
B cell: PAX5, MS4A1, CD19 …
Mac/mono markers: ADGRE1(또는 EMR1 = F4/80), LYZ (mouse: Lyz1), CD68, ITGAM (CD11b)…
Monocyte: CD14
NK: GZMA, GZMB, NKG7, GNLY, KLRD1…
Neutorphil: CSF3R, TREM1, S100A8
DC: ITGAX(CD11C), FLT3, CD1C, CLEC10A…
Platelet: PPBP, PF4
'''

# Let's draw each cell marker for FeatureFlot
FeaturePlot(pbmc_tutorial, c("CD3E", "MS4A1","LYZ", "EMR1", "CD14", "GZMA", "TREM1", "FLT3", "PPBP"), dims = c(1, 2))
# I assorted the locate according to the cluster's number on my own(Refer to Rplot05-2.jpg)
# writer: 0, 2, 4, 6 = T / 3 = B / 1, 5, 7 = Myeloid 계열 / 1번은 Monocyte / 6번은 NK / 7번은 DC / 8번은 Platelet으로 보입니다.
# re: where is the neurophil?? --;;

# AverageExpression, FindAllMarkers; to identify the more detailed features of clstrs
# and then... make a EXCEL file (no!!!!!!!!!!!!!!!!!!!!

Ave <- AverageExpression(pbmc_tutorial) # the latest ver. can be used
write.csv(Ave, "Cluster_avr.csv")
# create an object 'Ave' and make this to the CSV file

# make an excel sheet after searching marker genes on clusters via the function 'FindAllMarkers'
# install BiocManager & limma (skip)
Marker <- FindAllMarkers(pbmc_tutorial) # it will calculate clusters, No. zero to eight
write.csv(Marker, "Cluster_marker.csv")
# 실수로 위에 있는 파일 생성할 때 Cluster_avr.csv랑 똑같은 이름 만들어서 덮어쓰기로 저장함. 141번 행부터 재실행했는데 이게 맞는 대처 행위인지는 모르겠다... 이 때 이후로 시트 파일에 유전자 명이 이상하게 나타나는 것 같음(일단 원본 자료랑 발현도 내림차순 시 상위 에 뜨는 유전자 이름이 달라.)

# cluster annotation through the spreadsheet
# once selecting(filtering) the genes you want, put top 7 of amount of exp in the blank, Search for a gene (the link is below. click the 'search' and put ih in.)
# https://panglaodb.se/index.html
# then we can see Barplot of cell clusters (Y-axis) and cell types (X-axis) where the gene is expressed !!! but don't believe too much !!!
