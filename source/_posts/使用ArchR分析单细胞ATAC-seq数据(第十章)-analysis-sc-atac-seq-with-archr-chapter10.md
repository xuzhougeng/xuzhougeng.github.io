---
title: 使用ArchR分析单细胞ATAC-seq数据(第十章)
date: 2020-07-12 23:15:33.294
updated: 2020-07-12 23:17:27.211
url: /archives/analysis-sc-atac-seq-with-archr-chapter10
categories: R
tags: ATAC-seq | 单细胞
---

# 第10章: 使用ArchR识别Peak

识别peak(Peak calling)是ATAC-seq数据分析的常规操作之一。因为每个细胞的scATAC-seq数据都只有两种状态（开放或不开放），所以我们不能从单个细胞中识别peak。因此，我们在之前的章节中先对细胞进行了分组（细胞聚类）。此外，我们还构建了拟混池重复(pseduo-bulk replicates)，用于评估被识别的peak的可重复性。

## 10.1 迭代重叠Peak合并流程

我们先介绍[2018年的迭代重叠Peak合并流程](https://science.sciencemag.org/content/362/6413/eaav1898)，然后逐一介绍其他策略（这些策略都存在着一些问题）。

### 10.1.1 等宽和不定宽peak

我们之所以选择501-bp的等宽peak，一方面是下游处理不再需要根据peak宽度进行标准化从而简化了计算，另一方面是ATAC-seq的大部分peak宽度都低于501-bp。如果使用不定宽peak，后续不同样本的peak合并就会变得特别复杂。而且，我们也没有觉得使用不定宽peak带来的潜在好处值得我们花费那么多精力。更进一步的说，无论你选择哪种，大部分的分析结果也都不会受到影响。

下面，我们会以几个含有不同peak的细胞类型为例，用来讲解不同peak合并方法的差异。

### 10.1.2: Raw Peak Overlap

**raw peak overlap**就是分析不同细胞类型的peak是否有重叠，如果有将其合并成单个较大的peak。在如下的图解中，我们会发现尽管细胞类型A和C的前两个peak没有直接相连，但是由于细胞类型B中第一个peak和细胞类型A和C的前两个peak相互重叠，最终使得三个peak合并成一个更大peak。除了这个问题外，如果你想跟踪peak的顶点，你要么记录每一个合并后的新peak的顶点，要么记录每个合并后peak对应的原peak的顶点。

最后，这种类型的peak合并方法通常用`bedtools merge`命令实现。

<img src="/upload/2020/07/peakCalling_RawOverlap-a52800f5365a473a9626ccd63aca0755.png" alt="img" style="zoom: 25%;" />

### 10.1.3: Clustered Overlap

**clustered overlap**方法将peak进行聚类，然后根据规则挑选其中一个胜者。我们可以使用`bedtools cluster`命令对peak进行聚类，然后挑选每个peak聚类中最显著的peak。根据我们的经验，它会漏掉附近较小的peak，导致总体peak数减少。

<img src="/upload/2020/07/peakCalling_ClusteredOverlap-fb6ee0d46e45431785f5b5a11ef1ebaf.png" alt="img" style="zoom:25%;" />

### 10.1.4 Iterative Overlap

**iterative overlap**尽量避免了以上提到的问题。首先根据显著性对peak进行排序。接着我们保留其中最显著的peak，并移除掉所有于最显著peak发生重叠的peak。然后，在所有余下的peak中，我们会重复上述步骤，直到没有peak为止。这个能够避免10.1.2提到的问题，同时还能使用固定宽度peak。

<img src="/upload/2020/07/peakCalling_IterativeOverlap-e01c21c57f4a4493adef39a398beb44b.png" style="zoom:25%;" />

### 10.1.5 不同方法对比

我们很容易从最终的peak集中看到以上这些方法的区别。我们认为，**iterative overlap**策略能够较少的代价得到最好的peak集。

<img src="/upload/2020/07/peakCalling_Comparison-c8b214a7b733438298af50c3ff7dca3c.png" alt="comparison" style="zoom:25%;" />

### 10.1.6 ArchR的工作方式

**iterative overlap**这种以分级形式对数据进行合并，尽可能保留所有细胞类型特异的peak。

想象一下，你有3种细胞类型，A，B和C，每一种细胞类型多有3个拟混池重复。ArchR使用函数`addReproduciblePeakSet()`实现**iterative overlap**的peak合并流程。首先，ArchR单独为每个拟混池重复鉴定peak。然后ArchR同时分析每个细胞类型的所有拟混池重复，进行第一次的**iterative overlap**过滤。需要注意的是，ArchR使用标准化的peak显著性分值矩阵去比较不同样本的peak的显著性。因为有报道称MACS2的显著性和测序深度成比例关系，所以无法直接比较不同样本的peak显著性。在第一轮的**iterative overlap**过滤后，ArchR检查所有拟混池重复中每个peak的可重复性，只保留符合参数`reproducibility`的peak。最后，3种细胞类型(A, B, C)每一个都有单独的合并后的peak集。

接着，我们重复以上步骤，用于合并细胞类型A, B 和C的peak集。在这一步，我们还会根据不同的细胞类型对peak的显著性进行标准化，然后执行**iterative overlap**过滤。

最终我们会得到单个固定宽度的peak数据集。

### 10.1.7 如何更改策略

加入你不想用**iterative overlap**合并策略，那么你可以跳过`addReproduciblePeakSet()`，或者使用`ArchRProj <- addPeakSet()`添加自己的peak数据集

## 10.2 使用w/MACS2鉴定peak

如上所述，我们使用ArchR的`addReproduciblePeakSet()`得到可重复的peak。默认情况下，ArchR会使用MACS2鉴定peak，此外ArchR也实现它原生的peak鉴定工具，用在MACS2无法安装的情况（例如，我们无法在Windows上安装MACS2），该可选peak鉴定方法会在下一节介绍。

为了使用MACS2鉴定peak，ArchR必须能够找到MACS2的执行路径。ArchR会先在你的`PATH`环境变量进行查找。如果找不到，ArchR会尝试确认你是否用`pip`或`pip3`安装过MACS2. 如果以上都没有陈宫，ArchR就会放弃并输出报错信息。如果你安装了MACS2，但是ArchR找不到它，那么你需要在`addReproduciblePeakSet()`函数中设置`pathToMacs2`参数。

```r
pathToMacs2 <- findMacs2()
```

这是成功的信息

```bash
## Searching For MACS2..
## Found with $path!
```

如下是失败信息

```bash
Searching For MACS2..
Not Found in $PATH
Not Found with pip
Not Found with pip3
Error in findMacs2() : 
  Could Not Find Macs2! Please install w/ pip, add to your $PATH, or just supply the macs2 path directly and avoid this function!
```

当你有MACS2的位置信息后，我们接着就可以构建可重复的合并peak数据集(约5到10分钟)。为了避免细胞过少的拟混池重复的影响，我们可以通过设置`peaksPerCell`参数来限制每个细胞的peak数上限。这能避免细胞数过少的peak给最终合并的peak集贡献低质量的peak。此外，`addReproduciblePeakSet()`还有许多其他的参数可以调整，用`?addReproduciblePeakSet`可以看到更多信息。

每个`ArchRProject`项目只能包括一个peak集。我们最后将`addReproduciblePeakSet()`的输出赋值给我们想要的`ArchRProject()`。如果你想要测试不同的peak集，你必须先保存一份`ArchRProject()`并拷贝Arrow文件。尽管这会消耗更多的硬盘存储，但这是Arrow文件的会保存peak矩阵信息，因此无法避免。

```r
projHeme4 <- addReproduciblePeakSet(
    ArchRProj = projHeme4, 
    groupBy = "Clusters2", 
    pathToMacs2 = pathToMacs2
)
```

我们使用`getPeakSet()`函数以`GRanges`对象提取peak集。该peak集对象包括每个peak最初来源信息。然而这并不是意味着该peak只是在那个组中出现，而是意味着那组中的该peak拥有最高的标准化后显著性。

```r
getPeakSet(projHeme4)
```

## 10.3鉴定Peaks w/TileMatrix

如上所述，ArchR提供了原生的peak鉴定工具。尽管经过测试，该工具效果和MACS2表现差不多，但除非实在是没有办法，否则都不建议使用这个原生方法。

ArchR的原生peak鉴定工具基于500-bp的`TileMatrix`鉴定peak，可以通过设置`addReproduciblePeakSet()`函数的`peakMethod`参数来进行调用。注意，我们没有将输出保存在`projHeme4`对象中，因为我们这里只是测试而已，所以不打算保留这个peak集。将其保存在`ArchRProject`对象会覆盖之前已经存放在`projHeme4`的peak集。

```r
projHemeTmp <- addReproduciblePeakSet(
    ArchRProj = projHeme4, 
    groupBy = "Clusters2",
    peakMethod = "Tiles",
    method = "p"
)
```

我们同样可以检查这个合并后的peak集

```r
getPeakSet(projHemeTmp)
```

### 10.3.1 比较两种方法

我们可以根据重叠peak的百分比等指标来比较MACS2和ArchR自带的`TileMatrix`方法的差异。

1，我们检查MACS2鉴定的peak和`TileMatrix`鉴定的peak的重叠情况

```r
length(subsetByOverlaps(getPeakSet(projHeme4), getPeakSet(projHemeTmp))) / length(getPeakSet(projHeme4))
### 0.9627246
```

2，我们检查`TileMatrix`鉴定的peak和MACS2鉴定的peak的重叠情况。我们发现这个重叠比例并不高。

```r
length(subsetByOverlaps(getPeakSet(projHemeTmp), getPeakSet(projHeme4))) / length(getPeakSet(projHemeTmp))
### 0.7533365
```

如果我们增加peak的边界(从500-bp提高到1000)，MACS2鉴定的peak和`TileMatrix`鉴定的peak的重叠比例几乎没变化

```r
length(subsetByOverlaps(resize(getPeakSet(projHeme4), 1000, "center"), getPeakSet(projHemeTmp))) / length(getPeakSet(projHeme4))
### 0.9676687
```

但是增加了`TileMatrix`的peak和MACS2找到的peak的重叠比例。

```r
length(subsetByOverlaps(getPeakSet(projHemeTmp), resize(getPeakSet(projHeme4), 1000, "center"))) / length(getPeakSet(projHemeTmp))
### 0.8287639
```

来自译者的话: 我们从中可以得出，`TileMatrix`能鉴定的peak，MACS2也差不多都能鉴定到，但是MACS2鉴定出来的peak，未必`TileMatrix`也能鉴定。此外，`TileMatrix`找到的一些peak和MACS2找到的peak也比较近，因此通过延长peak的宽度，就能提高比例。

> 在10.1这个这个章节中，作者用到了菊花链拓扑(daisy chain)这个名词用来描述不同peak相互连接的情况。考虑到这个菊花链并非是一个常用名词，因此我通过意译来翻译这些名词。

## 10.4 增加Peak矩阵

我们可以用`saveArchRProject()`函数保存原来的`projHeme4`对象。该`ArchRProject`包含来源于MACS2的合并后peak集。

```r
saveArchRProject(ArchRProj = projHeme4, outputDirectory = "Save-ProjHeme4", load = FALSE)
```

为了方便下游分析，我们可以创建新的`ArchRProject`命名为`projHeme5`，并增加新的矩阵，该矩阵记录这新的合并后peak集的insertion计数。

```r
projHeme5 <- addPeakMatrix(projHeme4)
```

于是，`projHeme5`现在又多了一个新的矩阵，名为"Peakmatrix"，这是和"GeneScoreMatrix"和"TileMatrix"类似的**专用矩阵变量名**。正如之前所说，每个`ArchRProject`对象只能有一个peak集和一个`PeakMatrix`。当然，你可以不同命名构建无限个自定义的特征矩阵，但是`PeakMatrix`只能有一个，它就是用来保存来自于peak集的insertion计数矩阵。

```r
getAvailableMatrices(projHeme5)
## [1] "GeneIntegrationMatrix" "GeneScoreMatrix" "PeakMatrix"
## [4] "TileMatrix"
```






