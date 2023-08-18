---
title: motifmatchr:在R语言中分析peak中里是否有motif匹配
date: 2019-08-17 15:22:37.235
updated: 2019-09-10 22:57:26.142
url: /archives/motifmatchr-peak-motif-search-in-R
categories: 生信软件工具箱
tags: ATAC-seq | 表观组 | ChIP-seq
---


`motifmatchr`的作用就是分析众多的序列和众多的motifs, 从中找到哪个序列包含哪个motif. 它的核心函数就是`matchMotifs`，最大特点就是快，因为它用的是MOODS C++库用于motif匹配。

尽管Bioconductor上也有很多工具能够做motif匹配，比如说`Biostrings::mathcPWM`, `TFBSTools::searchSeq`，但是motifmatchr更适合于分析许多不同的序列包含许多不同的motif。例如，当分析ChIP-seq或者ATAC-seq数据时, 你可能会想知道在哪个peak里有哪种类型的motif.

## R包安装和加载

```r
library(motifmatchr)
library(GenomicRanges)
```

## matchMotifs

motifmatchr的核心函数是`matchMotifs`，因此了解这个函数的输入数据是什么格式就行了。必须的两个输入是

- 位置权重矩阵(position weight matrices, PWM)或位置频率矩阵(position frequency matrices, PFM), 保存在PWMatrix, PFMatrix, PWMatrixList或PFMatrixList
- 一组基因组范围(GenomicRanges或RangedSUmmarizedExperiment对象)或一组序列(DNAStringSet, DNAString 或 简单的字符串向量)

### PWM或PFM矩阵

Motif的PWM或PFM矩阵可以从[JASPAR](http://jaspar.genereg.net), [CIS-BP](http://cisbp.ccbr.utoronto.ca)下载。

例如，我从CIS-BP下载拟南芥Arabidopsis_thaliana_2019_08_17_2_32_am.zip，解压缩后里有一个`pwm_all_motifs`文件夹，里面的文本文件存放的就是我所需要的PWM矩阵,  下一步根据这些文件构建出matchMotifs的输入对象

```r
motif_dir <- "Arabidopsis_thaliana_2019_08_17_2_32_am/pwms_all_motifs"
PWList <- list()
for (file in list.files(motif_dir, pattern = ".txt")){
  df <- read.table(file.path(motif_dir, file), 
                   header = T,
                   row.names = 1)
  mt <- as.matrix(df)
  if (nrow(mt) ==0) next
  motif_id <- substr(file, 1,6)
  PWList[[motif_id]] <- PWMatrix(ID = motif_id, profileMatrix = t(mt))
}

PWMatrixLists <- do.call(PWMatrixList,PWList)
```

对于JASPAR，我们有更加方便的提取方法

```r
library(JASPAR2018)
species <- "Arabidopsis thaliana"
collection <- "CORE"
opts <- list()
opts["species"] <- species
opts["collection"] <- collection
out  <- TFBSTools::getMatrixSet(JASPAR2018, opts)

if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
  names(out) <- paste(names(out), TFBSTools::name(out), 
                      sep = "_")
motif <- out
```

JASPAR的motif比较可靠，因此motif相对比较少。

### 给定范围或序列

这部分的输入就比较容易获取，你可以提供MACS2的输出BED，利用`rtracklayer::import.bed()`读取构建成一个GRanges对象。因为是提供的GRanges对象，那么还需要额外设置一个参数`genome`, 利用`Biostrings::readDNAStringSet()`读取一个参考基因组序列就行了。

或者用bedtools先根据bed文件提取数据，然后利用`Biostrings::readDNAStringSet()`读取

### 示例数据

我们以包中提供的数据为例，进行演示

加载实例的motif和GRanges对象

```r
# load some example motifs
data(example_motifs, package = "motifmatchr") 

# Make a set of peaks
peaks <- GRanges(seqnames = c("chr1","chr2","chr2"),
                 ranges = IRanges(start = c(76585873,42772928,100183786),
                                  width = 500))

```

获取motif在peak中的位置

```r
motif_ix <- matchMotifs(example_motifs, peaks, genome = "hg19")
```

用`motifMatches`函数提取匹配矩阵

```r
motifMatches(motif_ix)
```

输出结果是一个稀疏矩阵

```r
3 x 3 sparse Matrix of class "lgCMatrix"
     MA0599.1_KLF5 MA0107.1_RELA MA0137.3_STAT1
[1,]             |             |              .
[2,]             |             .              |
[3,]             |             .              |
```

其中的`.`就表示存在motif匹配。

我们还可以提取motif在peak中的位置

```r
# Get motif positions within peaks for example motifs in peaks 
motif_ix <- matchMotifs(example_motifs, peaks, genome = "hg19",
                         out = "positions") 
```

### 其他参数

除了必须的motif信息和你的序列信息输入外，还有一些其他的参数可以做一些设置。

- 背景核苷酸频率, 默认是`bg=subject`, 也就是你的输入数据作为背景，也可设置成`genome`或`even`
- P值: 用于判断匹配是否足够可信的参数，默认是0.00005，基本上不用修改
- 输出信息: `matchMotifs`可以返回三种可能输出，matches, scores 和 positions

总的来说，这个R包是一个比较简单的工具，比较容易上手使用。