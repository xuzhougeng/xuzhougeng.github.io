---
title: snapATAC使用后的主观评测
date: 2021-12-30 06:05:11.034
updated: 2021-12-30 06:05:11.034
url: /archives/snapatac-subjective-evaluation
categories: 生信软件工具箱
tags: 单细胞
---

2019年，一篇发表在Genome Biology上题为「Assessment of computational methods for the analysis of single-cell ATAC-seq data」的文章对目前一些常见的单细胞ATAC-seq分析工具进行了评测，其中SnapATAC是当时文章认为比较好的工具

> Despite variation across methods and datasets, SnapATAC, Cusanovich2018, and cisTopic outperform other methods in separating cell populations of different coverages and noise levels in both synthetic and real datasets. Notably, SnapATAC is the only method able to analyze a large dataset (> 80,000 cells).

我曾经测试过这个工具，发现了其中的一个bug，后面就没有继续在用了。在2021年的最后几天，我又尝试使用了这个工具。如果你也在分析single cell ATAC-seq数据，想尝试snapATAC的话，希望这篇文章对你有所帮助。

1. snapATAC有一个配套的工具, snapTools, 它能够将fastq文件处理成所需的snap文件，或者将CellRanger-ATAC输出BAM/TSV文件转成snap文件。对于CellRanger-ATAC的数据转换，我不推荐使用BAM文件，因为他需要先给bam文件里read name添加barcode信息，然后按照read name进行排序，消耗的时间非常惊人。
2. snapATAC虽然能够处理大于80k的细胞，但是这和你选择的bin-size有关。snapATAC在分析的时候还是会将bin-by-cell的矩阵加载到内存中，因此内存还是一个比较重要的因素。我分析的时候，当bin-size比较小时，就会因为内存占用过多被kill。
3. snapATAC的可视化基于plot3D, 而不是ggplot2, 如果想要在原先图基础上进行修改，不如ggplot2那么方便。
4. snapATAC的聚类分析默认用的是R-igraph, leiden算法需要安装Python模块。R-igraph无法通过修改resolution来改变类群的数目，而leiden可以。不过，我好奇，为什么snapATAC不考虑使用Seurat的聚类算法？

