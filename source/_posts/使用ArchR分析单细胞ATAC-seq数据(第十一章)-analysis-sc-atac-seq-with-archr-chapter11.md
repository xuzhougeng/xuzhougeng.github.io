---
title: 使用ArchR分析单细胞ATAC-seq数据(第十一章)
date: 2020-07-19 06:58:51.204
updated: 2020-07-19 06:58:51.204
url: /archives/analysis-sc-atac-seq-with-archr-chapter11
categories: R
tags: ATAC-seq | 单细胞
---

# 第11章: 使用ArchR鉴定标记Peak

在讨论基因得分(gene score)这一章中，我们已经介绍了标记特征的鉴定。相同的函数`getMakerFeatures()`也能够用于从`ArchRProject`任意矩阵中鉴定标记特征。所谓的标记特征指的是相对于其他细胞分组唯一的特征。这些特征能帮助我们理解类群或者细胞类型特异的生物学现象。在这一章中，我们会讨论如何使用该功能鉴定标记peak。

## 11.1: 使用ArchR鉴定标记Peak

通常而言，我们想知道哪些peak是某个聚类或者某一些聚类所特有的。在ArchR中，这可以通过设置`addMarkFeatures()`函数的`useMatrix="PeakMatrix"`来实现（无需监督）。

首先，我们需要再看一眼`projHeme5`中有哪些细胞类型，以及它们的相对比例

```r
table(projHeme5$Cluster2)
```

现在，让我们调用`getMarkerFeatures`参数，并设置`useMatrix="PeakMatrix"`. 此外，为了降低不同细胞组之间的数据质量对结果的影响，我们可以设置`bias`参数，其中`bias = c("TSSEnrichment", "log10(nFrags)")`就是用来避免TSS富集和每个细胞的fragment数对结果的影响。

```r
markersPeaks <- getMarkerFeatures(
    ArchRProj = projHeme5, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
```

`getMarkerFeatures()`函数返回一个`SummarizedExperiment`对象，该对象包含一些不同的`assays`

```r 
markerPeaks
```

接着，我们可以用`getMarkers`函数从输出的`SummarizedExperiment`对象中提取我们感兴趣的部分。默认情况下，它会返回一个包含多个`DataFrame`的列表，不同的`DataFrame`表示来自不同的细胞分组。

```r
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
```

如果我们对特定的一个细胞分组感兴趣，我们可以用`$`进行提取。

```r
markerList$Erythroid
```

除了返回一个包含多个`DataFrame`的列表外，我们还可以用`getMarkers()`返回一个`GRangesList`，只要设置`returnGR=TRUE`即可。

```r
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markerList
```

这个`GRangesList`同样可以用`$`提取特定细胞组的结果，返回的是一个`GRanges`对象

```r
markerList$Erythroid
```

## 11.2 在ArchR中绘制Marker Peaks

ArchR提供了许多绘图函数用于`getMarkerfeatuers()`返回的`SummarizedExperiment`对象的可视化。

### 11.2.1 Marker Peak Heatmap

`markerHeatmap`能以热图的形式展示标记Marker Peak（或者其他`getMarkerFeatures()`输出的特征）

```r
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
```

使用`draw`函数绘制结果

```r
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)
```

![Marker Peak Heatmap](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/07/image-1459d659891f41b4a8f488ceb1fc8ee1.png)

`plotFDF()`函数能够以可编辑的矢量版本保存图片

```r
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)
```

### 11.2.2 Marker Peak MA和火山图

除了绘制热图，我们也可以为每个细胞分组绘制MA或者火山图(volcano)。这些图可以用`markerPlot()`函数绘制。对于MA图，需要设置参数`plotAs="MA"`. 我们以"Erythroid"细胞分组为例，设置参数`name = "Erythroid"`

```r
pma <- markerPlot(seMarker = markersPeaks, name = "Erythroid", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
pma
```

![MA图](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/07/image-904ba33c14ee4f6bab904eb1958a9732.png)

同样的，只要设置`plotAs="Volcano"`就可以绘制火山图

```r
pv <- markerPlot(seMarker = markersPeaks, name = "Erythroid", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
pv
```

![volcano plot](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/07/image-07e8fc062e764598a6ed6035b234baa5.png)

`plotFDF()`函数能够以可编辑的矢量版本保存图片。

```r
plotPDF(pma, pv, name = "Erythroid-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
```

### 11.2.3 Browser Tracks的Marker Peak

此外，我们在基因组浏览器上检查这些peak区域，只需要为`plotBrowserTrack()`函数的`features`参数传入 相应的peak区间。这会额外在我们的ArchR track图的下方以BED形式展示marker peak区域。

```r
p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "Clusters2", 
    geneSymbol = c("GATA1"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Erythroid"],
    upstream = 50000,
    downstream = 50000
)
```

我们使用`grid::grid.draw()`绘制结果

```r
grid::grid.newpage()
grid::grid.draw(p$GATA1)
```

![browser track](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/07/image-f86bc96dc4354f56aeba7db49c4fb4c9.png)

`plotFDF()`函数能够以可编辑的矢量版本保存图片。

```r
plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
```

## 11.3 组间配对检验

标记特征鉴定是一种特别的差异表达检验。此外，使用相同的`getMarkerFeatures()`函数也能实现标准化的差异分析。我们只需要设置`useGroup`为一组细胞，然后设置`bgdGroup`为另一组细胞即可。这就可以对给定两组进行差异分析。在这些差异分析中，在`useGroups`比较高的peak的倍数变化值是正数，在`bgdGroups`比较高的peak则是由负的倍数变化值。

这里，我们对"Erythroid"与"Progenitor"细胞组进行配对检验。

```r
markerTest <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Erythroid",
  bgdGroups = "Progenitor"
)
```

使用`markerPlot()`函数可以绘制MA或者火山图。MA图需要设置`plotAs="MA"`

```r
pma <- markerPlot(seMarker = markerTest, name = "Erythroid", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma
```

![group marker peak MA plot](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/07/image-e9db62e1bdbc4c709d979c08dbef8cb6.png)

火山图需要设置`plotAs="Volcano"`

```r
pv <- markerPlot(seMarker = markerTest, name = "Erythroid", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
pv
```

![group marker peak volcano](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/07/image-e4610c5ece3c47718f4514cdc17cd85b.png)

`plotFDF()`函数能够以可编辑的矢量版本保存图片。

```r
plotPDF(pma, pv, name = "Erythroid-vs-Progenitor-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
```

在后续章节中，我们还会进行差异分析，因为会在我们的差异开放的peak中搜索富集的motif。