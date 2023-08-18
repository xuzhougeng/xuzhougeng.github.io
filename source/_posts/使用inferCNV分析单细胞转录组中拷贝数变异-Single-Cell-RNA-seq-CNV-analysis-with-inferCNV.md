---
title: 使用inferCNV分析单细胞转录组中拷贝数变异
date: 2019-11-18 16:25:09.924
updated: 2019-11-18 16:25:09.924
url: /archives/Single-Cell-RNA-seq-CNV-analysis-with-inferCNV
categories: 生信软件工具箱
tags: 转录组 | 单细胞
---



`inferCNV`用与探索肿瘤单细胞RNA-seq数据，分析其中的体细胞大规模染色体拷贝数变化(copy number  alterations, CNA), 例如整条染色体或大片段染色体的增加或丢失(gain or deletions)。工作原理是，以一组"正常"细胞作为参考，分析肿瘤基因组上各个位置的基因表达量强度变化. 通过热图的形式展示每条染色体上的基因相对表达量，相对于正常细胞，肿瘤基因组总会过表达或者低表达。

`inferCNV`提供了一些过滤参数，通过调整参数来降低噪音，更好的揭示支持CNA的信号。此外`inferCNV`还包括预测CNA区间的方法以及根据异质性模式定义细胞类群的方法。

## 软件安装

尽管`inferCNV`是一个R包，但是在安装inferCNV之前还需要先下载安装[JAGS](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/) ，好在它有Windows，MacOS和Linux版本，所以`inferCNV`在各个平台都能用。

Windows和MacOS的JAGS容易安装，而Linux的JAGS需要编译

```bash
# 手动安装BLAS和LAPACK不推荐
# yum install blas-devel lapack-devel
tar xf JAGS-4.3.0.tar.gz 
cd JAGS-4.3.0
./configure --libdir=/usr/local/lib64
make -j 20 && make install 
```

安装R包

```R
install.packages("rjags")
if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
BiocManager::install("infercnv")
```

测试安装

```r
library(infercnv)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=system.file("extdata", "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv"),
                                    annotations_file=system.file("extdata", "oligodendroglioma_annotations_downsampled.txt", package = "infercnv"),
                                    delim="\t",
                                    gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=tempfile(), 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)
```

如果没有报错，就说明安装成功。

## 软件使用

### 准备输入文件

需要准备3个输入数据

1. 单细胞RNA-seq表达量的原始矩阵
1. 注释文件，记录肿瘤和正常细胞
1. 基因或染色体位置文件

第一个是Genes x Cells的表达矩阵(matrix)，行名是基因，列名是细胞编号。

| MGH54_P16_F12 | MGH54_P12_C10 | MGH54_P11_C11 | MGH54_P15_D06 | MGH54_P16_A03 | ...  |      |
| ------------- | ------------- | ------------- | ------------- | ------------- | ---- | ---- |
| A2M           | 0             | 0             | 0             | 0             | 0    | ...  |
| A4GALT        | 0             | 0             | 0             | 0             | 0    | ...  |
| AAAS          | 0             | 37            | 30            | 21            | 0    | ...  |
| AACS          | 0             | 0             | 0             | 0             | 2    | ...  |
| AADAT         | 0             | 0             | 0             | 0             | 0    | ...  |
| ...           | ...           | ...           | ...           | ...           | ...  | ...  |

第二个是样本注释信息文件，命名为"cellAnnotations.txt"。一共两列，第一列是对应第一个文件的列名，第二列是细胞的分组

```r
MGH54_P2_C12    Microglia/Macrophage
MGH36_P6_F03    Microglia/Macrophage
MGH54_P16_F12   Oligodendrocytes (non-malignant)
MGH54_P12_C10   Oligodendrocytes (non-malignant)
MGH36_P1_B02    malignant_MGH36
MGH36_P1_H10    malignant_MGH36
```

第三个是基因位置信息文件，命名为"geneOrderingFile.txt"。一共四列，第一列对应第一个文件的行名，其余三列则是基因的位置。**注**：基因名不能有重复

```r
WASH7P  chr1    14363   29806
LINC00115       chr1    761586  762902
NOC2L   chr1    879584  894689
MIR200A chr1    1103243 1103332
SDF4    chr1    1152288 1167411
UBE2J2  chr1    1189289 1209265
```

### 两步法

最复杂的工作就是准备输入文件，而一旦上述三个文件已经创建完成，那么分析只要两步以及根据结果对参数进行调整。

第一步，根据上述的三个文件创建`inferCNV`对象

```r
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=matrix, # 可以直接提供矩阵对象
                                    annotations_file="cellAnnotations.txt",
                                    delim="\t",
                                    gene_order_file="gene_ordering_file.txt",
                                    ref_group_names=c("normal"))
```

这一步的一个关键参数是`ref_group_name`, 用于设置参考组。假如你并不知道哪个组是正常，哪个组不正常，那么设置为`ref_group_name=NULL`, 那么`inferCNV`会以全局平均值作为基线，这适用于有足够细胞存在差异的情况。此外，这里的`raw_count_matrix`是排除了低质量细胞的count矩阵。

第二步，运行标准的`inferCNV`流程。

```r
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir",  # 输出文件夹
                             cluster_by_groups=T,   # 聚类
                             denoise=T, #去噪
                             HMM=T) # 是否基于HMM预测CNV
```

关键参数是`cutoff`, 用于选择哪些基因会被用于分析（在所有细胞的平均表达量需要大于某个阈值）。这个需要根据具体的测序深度来算，官方教程建议10X设置为0.1，smart-seq设置为1。你可以先评估下不同阈值下的保留基因数，决定具体值。`cluster_by_groups`用于声明是否根据细胞注释文件的分组对肿瘤细胞进行分群。

最终会输出很多文件在`out_dir`的目录下，而实际有用的是下面几个

- infercnv.preliminary.png : 初步的inferCNV展示结果（未经去噪或HMM预测）
- infercnv.png : 最终inferCNV产生的去噪后的热图.
- infercnv.references.txt : 正常细胞矩阵.
- infercnv.observations.txt : 肿瘤细胞矩阵.
- infercnv.observation_groupings.txt : 肿瘤细胞聚类后的分组关系.
- infercnv.observations_dendrogram.txt : NEWICK格式，展示细胞间的层次关系.

### 参数说明

`Infercnv::run`的参数非常之多，总体上分为如下几类

- 基本设置
- 平滑参数
- 基于肿瘤亚群的下游分析(HMM或non-DE-masking)
- 根据 BayesNet P(Normal) 过滤低可信度HMM预测结果
- 肿瘤亚群细分
- 去噪参数
- 离群值修剪
- 其他选项
- 实验性参数（不稳定）
- 差异表达分析的实验性参数

你可以按照具体的需求修改不同步骤的参数，例如聚类默认`cluster_by_groups=FALSE`会根据`k_obs_groups`聚类成指定的组数，而层次聚类方法用于计算组间相似度的参数则是`hclust_method`.

此外，设置`HMM=TRUE` 的计算时间会长于`HMM=FALSE`,因此可以先设置`HMM=FALSE`快速查看结果。

在运行过程中它会显示每个步骤的信息，官方文档给出了示意图帮助理解。

![InferCNV_procedure](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/11/InferCNV_procedure-e29ea5dab7b042f58bdc0600cb7ff489.png)

### 提取信息

`inferCNV`会输出一个" map_metadata_from_infercnv .txt"文件用于记录每个细胞的元信息，所有信息都可以从该文件中进行提取。或者使用`infercnv::add_to_seurat`将信息直接增加到原来的seurat对象中。

## 参考资料

- 软件安装: https://github.com/broadinstitute/inferCNV/wiki/Installing-infercnv
- 文件定义:  https://github.com/broadinstitute/inferCNV/wiki/File-Definitions 
- 运行inferCNV: https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV 
- **引用**: inferCNV of the Trinity CTAT Project.  https://github.com/broadinstitute/inferCNV

关于`inferCNV`的算法原理在如下几篇文章中有说明

- Anoop P. Patel, Itay Tirosh, et al. Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma. Science. 2014 Jun 20: 1396-1401
- Tirosh I et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science. 2016 Apr 8;352(6282):189-96
- Tirosh I et al. Single-cell RNA-seq supports a developmental hierarchy in human oligodendroglioma. Nature. 2016 Nov 10;539(7628):309-313. PubMed PMID: 27806376; PubMed Central PMCID: PMC5465819.
- Venteicher AS, Tirosh I, et al. Decoupling genetics, lineages, and microenvironment in IDH-mutant gliomas by single-cell RNA-seq. Science. 2017 Mar 31;355(6332).PubMed PMID: 28360267; PubMed Central PMCID: PMC5519096.
- Puram SV, Tirosh I, et al. Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor Ecosystems in Head and Neck Cancer. Cell. 2017 Dec 14;171(7):1611-1624.e24. PubMed PMID: 29198524; PubMed Central PMCID: PMC5878932.