---
title: 用rtracklayer读取和输出BigWig
date: 2019-08-13 19:58:24.434
updated: 2019-09-07 16:22:39.426
url: /archives/rtracklayer-BigWig
categories: R
tags: BigWig
---



BigWig是一个能用于加载到基因组浏览器上展示的格式。它的格式比较复杂，不适合直接阅读，通常由BedGraph文件转换而来。

在R语言中可以通过`rtracklayer`的`export.bw`输入和输出BigWig文件。

默认情况下，我们导入的数据集是GRanges格式

```r
test_path <- system.file("tests", package = "rtracklayer")
test_bw <- file.path(test_path, "test.bw")
gr <- import(test_bw)
gr 
```

输出信息如下

```r
gr
GRanges object with 9 ranges and 1 metadata column:
      seqnames    ranges strand |     score
         <Rle> <IRanges>  <Rle> | <numeric>
  [1]     chr2     1-300      * |        -1
  [2]     chr2   301-600      * |     -0.75
  [3]     chr2   601-900      * |      -0.5
  [4]     chr2  901-1200      * |     -0.25
  [5]     chr2 1201-1500      * |         0
  [6]    chr19 1501-1800      * |      0.25
  [7]    chr19 1801-2100      * |       0.5
  [8]    chr19 2101-2400      * |      0.75
  [9]    chr19 2401-2700      * |         1
  -------
  seqinfo: 2 sequences from an unspecified genome
```

我们可以通过设定`import`里的`as`参数，使其读取为RleList

```r
rle <- import(test_bw, as = "RleList")
rle
```

输出信息如下

```r
RleList of length 2
$chr2
numeric-Rle of length 243199373 with 5 runs
  Lengths:       300       300       300       300 243198173
  Values :        -1     -0.75      -0.5     -0.25         0

$chr19
numeric-Rle of length 59128983 with 6 runs
  Lengths:     1500      300      300      300      300 59126283
  Values :        0     0.25      0.5     0.75        1        0
```

如果要输出BigWig文件，也只要保证提供的对象是GRanges或RleList， 只不过要保证一定要有score列表示信号强度。

```r
export.bw(rle, "test.bw")
```

假如要用R语言输出下面这个GRanges对象为BigWig

```r
gr <- GRanges(
        seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges=IRanges(1:10, end=10),
        strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        seqlengths=c(chr1=11, chr2=12, chr3=13))
```

那么思路就是

1. 获取染色体的长度
1. 将染色体按照一定宽度分窗, 可用函数为`tile`和`slidingWindows`
1. 统计每个窗口的read数, 需要用到函数`coverage`, `Views`和`ViewSums`
1. 输出BigWig

代码实现:

```r
# 获取染色体长度
chromSizes <- GRanges(c("chr1","chr2","chr3"), IRanges(1,c(11,12,13)))
# 按照宽度2进行分窗
windows <- tile(chromSizes, width = 2)
# 统计每个窗口的read数
cvg <- coverage(gr)
cov_by_wnd <- Views(cvg, windows)
sum_by_wnd <- viewSums(cov_by_wnd)
for(i in seq_along(windows)){
  mcols(windows[[i]])['score'] <- sum_by_wnd[[i]]
}
# 输出
seqlengths(windows) <- c(11,12,13)
export.bw(unlist(windows), "tmp.bw")
```

