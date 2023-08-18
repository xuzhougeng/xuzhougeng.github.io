---
title: Seurat的scATAC-seq分析流程
date: 2019-08-31 10:10:38.344
updated: 2019-09-01 13:26:12.556
url: /archives/Analysis-scATAC-seq-with-scRNA-seq-in-Seurat
categories: 生信软件工具箱
tags: Seurat | 转录组 | 单细胞
---


Seurat 3.X版本能够整合scRNA-seq和scATAC-seq, 主要体现在：

- 基于scRNA-seq的聚类结果对scATAC-seq的细胞进行聚类
- scRNA-seq和scATAC-seq共嵌入(co-embed)分析

整合步骤包括如下步骤:

1. 从ATAC-seq中估计RNA-seq表达水平，即从ATAC-seq reads定量基因表达活跃度
1. 使用LSI学习ATAC-seq数据的内部结构
1. 鉴定ATAC-seq和RNA-seq数据集的锚点
1. 数据集间进行转移，包括聚类的标签，在ATAC-seq数据中推测RNA水平用于共嵌入分析

## 数据下载

测试数据下载地址:

- [scATAC-seq](http://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5): h5格式

- [scATAC-seq metadata](http://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv): csv文件
- [scRNA-seq](http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5): h5格式

可以复制下载链接到浏览器下载，也可以直接在R语言用`download.file`中进行下载。

```R
# 下载peak
atac_peak <- "http://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
atac_peak_file <- "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
download.file(atac_peak, atac_peak_file)

# 下载singlecell.csv
singlecell <- "http://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv"
singlecell_file <- "atac_v1_pbmc_10k_singlecell.csv"
download.file(singlecell, singlecell_file)

# 下载rna-seq
rna_seq <- "http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
rna_seq_file <- "pbmc_10k_v3_filtered_feature_bc_matrix.h5"
download.file(rna_seq, rna_seq_file)

# 下载GTF
gtf <- "ftp://ftp.ensembl.org/pub/grch37/release-82/gtf/homo_sapiens/Homo_sapiens.GRCh37.82.gtf.gz"
gtf_file <- "Homo_sapiens.GRCh37.82.gtf.gz"
download.file(gtf, gtf_file)
```

## 基因活跃度定量

首先，先将peak矩阵转成基因活跃度矩阵。Seurat做了一个简单的假设，基因活跃度可以通过简单的将落在基因区和其上游2kb的count相加得到基因活跃度，并且这个结果Cicero等工具返回gene-by-cell矩阵是类似的。

```R
library(Seurat)
library(ggplot2)
# 读取peak
peaks <- Read10X_h5(atac_peak_file)
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, 
                                            annotation.file = gtf_file, 
                                            seq.levels = c(1:22, "X", "Y"), 
                                            upstream = 2000, 
                                            verbose = TRUE)
```

activity.matrix是一个`dgCMatrix`对象，其中行为基因，列为细胞。因此如果对于Cicero的输出结果，只要提供相应的矩阵数据结构即可。


## 设置对象

下一步，我们要来设置`Seurat`对象，将原始的peak counts保存到assay中，命名为"ATAC"

```r
pbmc.atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "10x_ATAC")
```

增加基因活跃度矩阵，命名为"ACTIVITY".  

```r
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
```

增加细胞的元信息，该信息来自于scATAC-seq的CellRanger处理的singlecell.csv

```r
meta <- read.table(singlecell_file, sep = ",", header = TRUE, row.names = 1, 
    stringsAsFactors = FALSE)
meta <- meta[colnames(pbmc.atac), ]
pbmc.atac <- AddMetaData(pbmc.atac, metadata = meta)
```

过滤掉scATAC-seq数据中总count数低于5000的细胞，这个阈值需要根据具体实验设置

```R
pbmc.atac <- subset(pbmc.atac, subset = nCount_ATAC > 5000)
pbmc.atac$tech <- "atac"
```

## 数据预处理

为了找到scATAC-seq数据集和scRNA-seq数据集之间的锚定点，我们需要对基因活跃度矩阵进行预处理

设置pbmc.atac的默认Assay为"ACTIVITY"， 然后寻找高表达的基因，对基因活跃度矩阵进行标准化和Scale。

```R
DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- FindVariableFeatures(pbmc.atac)
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac)
```

同样的，我们还需要对peak矩阵进行处理。这里我们用的隐语义(Latent semantic indexing, LSI)方法对scATAC-seq数据进行降维。该步骤学习scRNA-seq数据的内部结构，并且在转换信息时对锚点恰当权重的决定很重要。

根据 [Cusanovich et al, Science 2015](https://science.sciencemag.org/content/348/6237/910/tab-pdf)提出的LSI方法，他们搞了一个`RunLSI`函数。LSI的计算方法为TF-IDF加SVD。

我们使用在所有细胞中至少有100个read的peak，然后降维到50。该参数的选择受之前的scATAC-seq研究启发，所以没有更改，当然你你可以把它改了。

```r
DefaultAssay(pbmc.atac) <- "ATAC"
VariableFeatures(pbmc.atac) <- names(which(Matrix::rowSums(pbmc.atac) > 100))
pbmc.atac <- RunLSI(pbmc.atac, n = 50, scale.max = NULL)
#pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 1:50)
pbmc.atac <- RunTSNE(pbmc.atac, reduction = "lsi", dims = 1:50)
```

**注**: 要将pbmc.atac的默认assay切换成"ATAC",  非线性降维可以选择UMAP或者t-SNE。

我们之前使用过Seurat对scRNA-seq数据进行预处理和聚类，下载地址为[Dropbox](https://www.dropbox.com/s/3f3p5nxrn5b3y4y/pbmc_10k_v3.rds?dl=1)。

```r
pbmc.rna <- readRDS("pbmc_10k_v3.rds")
pbmc.rna$tech <- "rna"
```

将scRNA-seq和scATAC-seq共同展示，对一些骨髓(myeloid)和淋巴(lymphoid)PBMC中比较常见的聚类，其实是能从图中看出来。

```r
p1 <- DimPlot(pbmc.atac, reduction = "tsne") + 
        NoLegend() + 
        ggtitle("scATAC-seq")
p2 <- DimPlot(pbmc.rna, group.by = "celltype", reduction = "tsne", label = TRUE, repel = TRUE) + 
        NoLegend() + 
        ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))
```

![标签转移前](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/without-transfer-b47eb2c142b047d58afcc64b5068115a.png)

现在，我们用`FindTransferAnchors`鉴定scATAC-seq数据集和scRNA-seq数据集的锚点，然后根据这些锚点将 10K scRNA-seq数据中鉴定的细胞类型标记转移到scATAC-seq细胞中。

```r
transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, 
                                        query = pbmc.atac, 
                                        features = VariableFeatures(object = pbmc.rna), 
                                        reference.assay = "RNA", 
                                        query.assay = "ACTIVITY", 
                                        reduction = "cca")
```

> Seurat比较的是scRNA-seq表达量矩阵和scATAC-seq中基因活跃度矩阵，利用CCA降维方法比较两者在scRNA-seq中的高变异基因的关系

为了转移细胞类群的编号，我们需要一组之前注释过的细胞类型，作为`TransferData`的 refdata 参数输入。`TransferData`本质上采用的是KNN算法，利用距离未知类型细胞的最近N个已知细胞所属的类型来定义该细胞。`weight.reduction`参数是用来选择设置权重的降维。

```r
celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = pbmc.rna$celltype, 
                                     weight.reduction = pbmc.atac[["lsi"]])
pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)
```

我们可以检查每个细胞预测得分的分布情况，选择性的过滤哪些得分比较低的细胞。我们发现超过95%的细胞得分大于或等于0.5.

```r
hist(pbmc.atac$prediction.score.max)
abline(v = 0.5, col = "red")
```

![得分分布](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/1567213591418-23ede2d0b7324a30b2aa10dda80ba7f8.png)

```r
table(pbmc.atac$prediction.score.max > 0.5)
# FALSE  TRUE 
#   383  7483 
```

之后，我们就可以在UMAP上检查预测的细胞类型的分布，检查细胞类型在scRNA-seq和scATAC-seq中的相对位置

```r
pbmc.atac.filtered <- subset(pbmc.atac, 
                             subset = prediction.score.max > 0.5)
pbmc.atac.filtered$predicted.id <- factor(pbmc.atac.filtered$predicted.id, 
                                          levels = levels(pbmc.rna))  # to make the colors match
p1 <- DimPlot(pbmc.atac.filtered, 
              group.by = "predicted.id",
              label = TRUE, 
              repel = TRUE) + 
        ggtitle("scATAC-seq cells") + 
        NoLegend() + 
        scale_colour_hue(drop = FALSE)
p2 <- DimPlot(pbmc.rna, group.by = "celltype", label = TRUE, repel = TRUE) +
        ggtitle("scRNA-seq cells") + 
        NoLegend()
CombinePlots(plots = list(p1, p2))
```

![预测后结果](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/1567213664297-0e4113573a6d47d081e5699e064a1624.png)

在转移细胞类型标记之后，我们就可以进行细胞特意水平上的下有分析。举个例子，我们可以去找一些某些细胞类型特异的增强子，寻找富集的motif。目前这些分析Seurat还不直接支持，还在调试中。

## 共嵌入(co-embedding)

最后，如果你想将所有的细胞一同展示，可以将scRNA-seq和scATAC-seq数据嵌入到相同的低维空间。

我们使用之前的锚点从scATAC-seq细胞中推断RNA-seq的值，后续分析就相当于两个单细胞数据的分析流程。

**注意:** 这一步只是为了可视化，其实不做也行。

选择用于估计的基因，可以是高变动基因，也可以是所有基因。

```r
# 只选择高变动的基因作为参考
genes.use <- VariableFeatures(pbmc.rna)
refdata <- GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use, ]
```

之后利用`TransferData`推断scATAC-seq在这些基因中的可能值，输出结果就是ATAC细胞对应的scRNA-seq矩阵

```r
imputation <- TransferData(anchorset = transfer.anchors, 
                           refdata = refdata, 
                           weight.reduction = pbmc.atac[["lsi"]])
# this line adds the imputed data matrix to the pbmc.atac object
pbmc.atac[["RNA"]] <- imputation
```

合并两个的结果，然后就是scRNA-seq的常规分析。

```r
coembed <- merge(x = pbmc.rna, y = pbmc.atac)
# 标准化分析
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
#coembed <- RunUMAP(coembed, dims = 1:30)
coembed <- RunTSNE(coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)
```

在t-SNE上绘制结果

```r
p1 <- DimPlot(coembed, reduction="tsne", group.by = "tech")
p2 <- DimPlot(coembed, reduction="tsne", group.by = "celltype", label = TRUE, repel = TRUE)
CombinePlots(list(p1, p2))
```

![tSNE plot](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/1567215158904-47bba1920748496bbaaf95e7a5ec773f.png)

从上面的结果中，你可能会发现某些细胞只有在一类技术中存在。举个例子，从巨噬细胞(megakaryocytes)到成熟的血小板细胞(pletelet)由于没有细胞核，无法被scATAC-seq技术检测到。我们可以单独展示每个技术，进行检查

```r
DimPlot(coembed, reduction="tsne", split.by = "tech", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()
```

![分别展示](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/1567215198865-9b93b77c817b445eb0d00bf1b0455821.png)

此外，你还可以发现B细胞前体类型附近有一类细胞只由scATAC-seq的细胞构成。通过展示这些细胞在CellRanger分析结果提供的黑名单位置所占read数，可以推测出这类细胞可能是死细胞，或者是其他scRNA-seq无法检测的人为因素。

```r
coembed$blacklist_region_fragments[is.na(coembed$blacklist_region_fragments)] <- 0
FeaturePlot(coembed, features = "blacklist_region_fragments", max.cutoff = 500)
```

![黑名单区](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/1567217384234-ef3b56c946154ef993515b57223f2c6b.png)

---

Seurat这种基于基因活跃度得分进行细胞类型预测，是否靠谱，开发者提供了如下几个证据

- 总体预测得分(`pbmc.atac$prediction.score.max`)高，意味着用scRNA-seq定义细胞类型比较可靠
- 我们可以在scATC-seq降维结果中
- 利用相同锚点的贡嵌入分析，发现两类形态能很好的混合
- 将ATAC-seq数据根据聚类结果构建pseduo bulk, 发现和真实的bulk数据近似

## 参考资料：

- <https://satijalab.org/seurat/v3.0/atacseq_integration_vignette.html>