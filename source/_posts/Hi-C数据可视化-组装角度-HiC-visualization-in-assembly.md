---
title:  Hi-C数据可视化-组装角度
date: 2019-12-08 22:49:33.209
updated: 2019-12-08 22:49:33.209
url: /archives/HiC-visualization-in-assembly
categories: 生信软件工具箱
tags: 组装 | Hi-C
---

# Hi-C数据可视化-组装角度

> 这里讨论HiC的可视化是从组装角度出发，也就是如何展示contig和contig之间的关系

Hi-C数据可视化（我所知）有下面几个

- [Juicerbox](https://github.com/aidenlab/Juicebox/wiki/Download): 可视化图形工具，需要.hic作为输入。
- [HiTC](http://bioconductor.org/packages/release/bioc/html/HiTC.html): R包工具
- [HiCPlotter](https://github.com/kcakdemir/HiCPlotter): Python工具
- [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/index.html): Python工具
- [HiGlass](higlass.io): 网页浏览器

无论是何种分析工具，最核心的步骤就是提供它们需要的输入数据，也就是将fastq转成所需要的输入形式。虽然本质上都是对contact matrix进行展示，但问题在于部分工具定义了其特有的存放格式，需要进行一些格式转换。

基于上述考虑，我在数据预处理上选择使用[HiC-Pro](/archives/HiC-Pro-An-optimized-and-flexible-pipeline-for-Hi-C-data-processing)，因为`HiC-Pro`本身就能够将其结果转成fithic,juicebox的输入，此外`HiCExplorer`提供了`hicCovertFormat`，能够将Hic-Pro转成cool/mcool/h5/homer/ginteractions等格式。

后续假定HiC-Pro输出的文件为三个文件

- input_20000.matrix: 原始上三角矩阵
- input_20000_abs.bed: 坐标信息
- input_20000_iced.matrix: ICE normalization的矩阵

## HiTC

HiTC是R包，能够直接读取HiC-pro的输出结果，对应的函数是`importC`

```R
# BiocManager::install("HiTC")
library(HiTC)
## Load Hi-C data
x <- importC("input_20000.matrix", xgi.bed = "input_20000_abs.bed")
show(x)
```

由于读取这一步时，它会创建NXN个 Map，如果有1,000个contig， 那么最终会得到1,000,000个map，这是一个非常可怕的计算量，R语言会因此而崩溃。因此对于contig级别的基因组而言，不应该用这个工具进行展示。

在可视化方面方面，主要是`mapC`，它接受一个HTClist作为输入

```r
sset <- reduce(x, chr=c("chr5","chr6","chr7"))
x90_500 <- HTClist(mclapply(sset, binningC,
binsize=500000, bin.adjust=FALSE, method="sum", step=1))
mapC(imr90_500)
```

考虑到矩阵以对角线对称，实际上只需要展示上三角即可

```r
hox <- extractRegion(x$chr6chr6, chr="chr6", from=50e6, to=58e6)
plot(hox, maxrange=50, col.pos=c("white", "orange", "red", "black"))
```

更多参考: <https://bioconductor.org/packages/release/bioc/vignettes/HiTC/inst/doc/HiC_analysis.pdf>


## HiCPlotter

HiCPlotter依赖Python2.7和对应的Numpy, Scipy,Matplotlib, 由于两年没有更新，因此GitHub里提到的版本都比较古老，经测试，Numpy=1.15.1，Scipy=1.1.0， Matplotlib=2.2.3 能够正常运行。

它支持HiC-pro的输出，无需做格式转换。

举个例子

```bash
python HiCPlotter.py \
    -f input_20000_iced.matrix \
    -o Exemple \
    -wg 1
    -r 20000 \
    -tri 1 \
    -bed input_20000_abs.bed \
    -n Test -chr chrX -ptr 1
# -r: 分辨率
# -f 输入矩阵
# -tri: HiC-Pro输入要设置为1
# -bed: HiC-Pro输入要指定bed
# -wg: 全基因组的交互
# -ptr: 绘制三角还是全部
```

HiCPlotter只能单个染色体或者从起始染色体到指定染色体的全基因组互作，并不能绘制给定任意多个染色体的互作关系，功能比较简单。

## HiCExplorer可视化

HiCExplorer开发deeptools的团队推出的HiC处理工具，包括格式转换，HiC矩阵质控和可视化。在功能上和性能上，HiCExplorer都比之前提及的工具强，但是如果只用需要的功能，它并不会很复杂。

HiCExplorer要求contact matrix以h5格式存放，因此需要对HiC-Pro的结果进行转换

```bash
hicConvertFormat \
	-m input_20000_iced.matrix \
	-o input_20000 \
	--bedFileHicpro input_20000_abs.bed \
	--inputFormat hicpro --outputFormat h5
```

输出的是input_20000.h5文件，该h5文件便是后续分析的矩阵输入文件

使用`hicPlotMatrix`绘制热图展示HiC数据

```bash
# 绘制单条染色体
hicPlotMatrix --matrix input_20000.h5 -out chrX_chrY \
    --region chrX 
# chrX和chrY的交互
hicPlotMatrix --matrix input_20000.h5 -out chrX_chrY \
    --region chrX \
    --region2 chrY
# 展示指定的2个染色体
hicPlotMatrix --matrix input_20000.h5 -out chrXchrY \
    --chromosomeOrder chrX chrY
```

更多信息参考: <https://hicexplorer.readthedocs.io/en/latest/content/tools/hicPlotMatrix.html>



