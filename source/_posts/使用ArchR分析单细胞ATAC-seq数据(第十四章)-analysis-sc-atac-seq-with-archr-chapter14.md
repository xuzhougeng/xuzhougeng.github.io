---
title: 使用ArchR分析单细胞ATAC-seq数据(第十四章)
date: 2021-01-17 16:29:32.657
updated: 2021-01-17 16:31:04.137
url: /archives/analysis-sc-atac-seq-with-archr-chapter14
categories: 生信软件工具箱
tags: ATAC-seq | 单细胞
---

# 第十四章 ArchR的足迹分析

转录因子(Transcripts factor, TF)足迹分析使得我们能够预测特定位点中TF的精确结合位置。这是因为该位置被TF结合避免了转座酶的切割，而TF结合位点的邻近位置处于开放状态。


![image.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/01/image-9c35998c679f4262b1d3f7fe73533199.png)

理想情况下，TF足迹分析需要在单个位置上分析从而确定TF的准确结合位置。但实际上，这需要非常高的测序深度，甚至超过混池ATAC-seq或者scATAC-seq的所有数据。为了解决这个问题，我们可以把和待预测的TF结合相关的Tn5插入位置进行合并。例如，我们可以提取所有包含CTCF motif的peak，制作一个全基因组的CTCF的聚合TF足迹。

为了保证足迹的可靠性，我们需要确保能够可靠的预测出目标TF所对应的结合位点。ArchR使用自带的`addMotifAnnotations()`函数对peak区域进行搜索，寻找能够匹配的DNA序列。考虑到motif的简并性，无法保证每个motif都有足够的peak。添加到`ArchRProject`的motif注释以二值矩阵表示(0=无motif, 1=有motif)。一旦你有了这些motif注释，ArchR使用`getFootPrints()`函数分析足迹，它以一个`ArchRProject`对象和一个`GenomicRanges`对象(记录motif的位置)作为输入。可以使用`getPositions()`函数从`ArchRProject`中提取这些位置。之后足迹可以使用`plotFootprints()`函数可视化。

或许更重要的是，ArchR的足迹分析能够抵消已知的Tn5插入序列偏好性。ArchR使用一个hexmer位置频率矩阵和一个目标Tn5插入位置上的k-mer频率矩阵来实现该功能。

![foot printing Methods](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/01/image-1c49ea2ad32e4f7b933132461b8233d0.png)


最终，该流程输出考虑到Tn5插入偏好性的足迹图。

ArchR支持motif足迹分析和用户提供特征的足迹分析，在后续都会讨论。

## 14.1 motif足迹分析

由于教程用到的数据集比较小，因此利用该数据集得到足迹并不是那么的清晰。使用更大的数据会得到较小变异的足迹。

在分析足迹时，我们要先获取和motif相关的所有位置。这可以通过`getPositions`函数完成。该函数有一个可选参数， `name`，用于传入`peakAnnotation`对象中我们想要获取位置的变量名。如果`name=NULL`, 那么ArchR会使用`peakAnnotation`槽(slot)的第一个条目(entry)。在下面的例子中，我们没有指定`name`, ArchR使用的第一个条目为CIS-BP motifs.

```r
motifPositions <- getPositions(projHeme5)
```

这会创建一个`GRangesList`对象，每个TF motif以不同的`GRanges`对象进行区分。

```r
motifPositions
# GRangesList object of length 870:
# $TFAP2B_1
# GRanges object with 16773 ranges and 1 metadata column:
# seqnames ranges strand | score
# |
# [1] chr1 852468-852479 + | 8.17731199746359
# [2] chr1 873916-873927 + | 8.32673820065588
# [3] chr1 873916-873927 - | 8.32673820065588
# [4] chr1 896671-896682 + | 9.96223327271814
# [5] chr1 896671-896682 - | 8.92408377606486
# … … … … . …
# [16769] chrX 153991101-153991112 + | 8.39549159740639
# [16770] chrX 154299568-154299579 + | 8.90119825654299
# [16771] chrX 154664929-154664940 - | 8.16690864294221
# [16772] chrX 154807684-154807695 + | 9.57636587154549
# [16773] chrX 154807684-154807695 - | 10.6117355833828
# ——-
# seqinfo: 23 sequences from an unspecified genome; no seqlengths
#
# …
# <869 more elements>
```

我们提取部分感兴趣的TF motifs用于展示。在我们提取"EBF1"会附带"SREBF1", 因此我们需要用`%ni%`显式将其过滤。`%ni%`函数是R自带函数`%in%`的相反函数。

```r
motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs
# [1] “GATA1_383” “CEBPA_155” “EBF1_67” “IRF4_632” “TBX21_780” “PAX5_709”
```

为了准确找到TF足迹，我们需要大量的reads。因此，细胞需要进行分组生成拟混池ATAC-seq谱才能用于TF足迹分析。这些拟混池谱之前在peak鉴定时就已经保存为分组覆盖文件。 如果没有在`ArchRProject`添加分组覆盖信息，则运行如下命令

```r
projHeme5 <- addGroupCoverages(ArchRProj = projHeme5, groupBy = "Clusters2")
```

在计算分组覆盖度后，我们可以为之前`getFootprints()`挑选的一组标记motif计算足迹。即便ArchR已经优化了足迹分析流程，我们也建议先对一部分motif分析足迹，而不是直接分析所有motif。 我们通过`positions`参数来选择motif。

```r
seFoot <- getFootprints(
  ArchRProj = projHeme5, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters2"
)
```

当我们获取了这些足迹，我们可以使用`plotFootprints()`函数进行展示。该函数能够同时以多种方式对足迹进行标准化。下一节，我们会讨论标准化和实际的足迹图。

## 14.2 Tn5偏好的足迹标准化

使用ATAC-seq数据分析TF足迹的一大挑战就是Tn5转座酶的插入序列偏好性，这会导致TF足迹的错误分类。为了降低Tn5插入偏好性的影响，ArchR识别每个Tn5插入位置附近的k-mer序列(k由用户提供，默认是6). 

对于该项分析，ArchR为每个拟混池识别单碱基分辨率的Tn5插入位点，将这些1-bp位点调整为k-bp窗口（-k/2和+（k/2-1）bp），然后使用`Biostrings`包中的`oligonucleotidefrequency(w=k, simplify.as="collapse")`函数创建k-mer频率表。然后，ArchR使用与`BSgenome`相关的基因组文件，以相同的函数计算出全基因组范围预期的k-mers。

为了计算拟混池足迹的插入偏差，ArchR创建了一个k-mer频率矩阵，该矩阵表示为从motif中心到窗口+/-N bp（用户定义，默认为250 bp）的所有可能k-mer。然后，遍历每个motif位点，ArchR将定位的k-mer填充到k-mer频率矩阵中。然后在全基因组范围内计算每个motif位置。利用样本的k-mer频率表，ArchR可以通过将k-mer位置频率表乘以观察/期望 Tn5 k-mer频率来计算预期的Tn5插入。

以上所有这些发生在`plotFootprints()`函数中。

### 14.2.1 减去Tn5偏好

一个标准化方式就是从足迹信号中减去Tn5偏好。该标准化方法通过设置`plotFootprints()`的`normMethod = "Subtract"`实现

```r
plotFootprints(
  seFoot = seFoot,
  ArchRProj = projHeme5, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)
```

默认，这些图保存在`ArchRProjec`t的`outputDirectory`。如果你需要绘制所有motif, 可以将其返回为`ggplot2`对象，需要注意这个`ggplot`对象会非常大。下面是一个从motif足迹中减去Tn5偏好信号的结果

![Footprints-Subtract-Bias](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/01/image-9445c244d46f4c60b77b93eaf5209d14.png)

### 14.2.2 除以Tn5偏好

第二种标准化方法就是将足迹除以Tn5偏好信号。该标准化方法通过设置`plotFootprints()`的`normMethod = "Divide"`实现

```r
plotFootprints(
  seFoot = seFoot,
  ArchRProj = projHeme5, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)
```

下面是一个从motif足迹中除以Tn5偏好信号的结果

![Footprints-Divide-Bias](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/01/image-16a63ea23fac4d31abf71b7b09d9a091.png)

### 14.2.3 无Tn5偏好标准化的足迹

尽管我们高度推荐将足迹根据Tn5序列插入偏好性进行标准化，当然你可以通过设置`plotFootprints()`的`normMethod = "None"`来省去标准化。

```R
plotFootprints(
  seFoot = seFoot,
  ArchRProj = projHeme5, 
  normMethod = "None",
  plotName = "Footprints-No-Normalization",
  addDOC = FALSE,
  smoothWindow = 5
)
```


![Footprints-No-Normalization](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/01/image-7fcb07f110c747c4ac55c3a71bde7f96.png)


## 14.3 特征足迹

除了motif足迹分析，ArchR还允许用户分析任意定义的特征数据集。为了对功能进行说明，我们将会使用`plotFootprints()`函数创建TSS插入谱(在之前数据质量控制一节中引入)。一个TSS插入谱本质上就是特殊的足迹。

我们在之前小节讨论过，足迹会用到来源于拟混池重复的分组覆盖文件。我们在之前peak鉴定时创建过这些文件。如果你没有在`ArchRProject`加入分组覆盖信息, 那么需要运行如下代码

```r
projHeme5 <- addGroupCoverages(ArchRProj = projHeme5, groupBy = "Clusters2")
```

我们接着创建一个没有经过Tn5偏好性校正的TSS插入谱。和之前分析一个主要不同是，我们设置了`flank=2000`, 将TSS向前向后分别延伸2000 bp.

```r
seTSS <- getFootprints(
  ArchRProj = projHeme5, 
  positions = GRangesList(TSS = getTSS(projHeme5)), 
  groupBy = "Clusters2",
  flank = 2000
)
```

我们接着用`plotFootprints()`对每个细胞分组绘制TSS插入谱。

```r
plotFootprints(
  seFoot = seTSS,
  ArchRProj = projHeme5, 
  normMethod = "None",
  plotName = "TSS-No-Normalization",
  addDOC = FALSE,
  flank = 2000,
  flankNorm = 100
)
```

![TSS-No-Normalization](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/01/image-b23beabc98094afd8d59ba25139d06e6.png)

