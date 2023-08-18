---
title: 如何展示MUMMER的结果
date: 2019-10-10 23:06:10.385
updated: 2019-10-10 23:06:10.385
url: /archives/Visualization-of-MUMMER-Result
categories: R
tags: 可视化
---



在2年之前我写过一篇教程介绍MUMmer软件的使用方法，可以通过[如何使用MUMmer比对大片段序列](/archives/Using-MUMmer-to-align-genome)阅读。

MUMmer作为一个比对工具，它的主要功能就是寻找两个序列的相似之处，至于如何展示结果，并不是它的主要目标。这篇文章将会介绍如何基于MUMMER的输出结果进行可视化。

首先是下载数据，我们用细菌的基因组作为案例

```bash
wget http://mummer.sourceforge.net/examples/data/H_pylori26695_Eslice.fasta
wget http://mummer.sourceforge.net/examples/data/H_pyloriJ99_Eslice.fasta
```

然后使用`nucmer`进行比对

```bash
nucmer H_pylori26695_Eslice.fasta H_pyloriJ99_Eslice.fasta
```

只保留1对1的最佳联配

```bash
delta-filter -1 out.delta > filter.delta
```

输出容易展示的信息

```bash
show-coords -T -q -H filter.delta > coord.txt
```

上面这几个程序的参数需要根据具体需求进行修改，并不是固定用法

后续就是在R语言里面继续绘图。 先加载数据并根据数据格式命名列

```r
df <- read.table("~/coord.txt", sep = "\t")

colnames(df) <- c("ref_start", "ref_end", "qry_start", "qry_end", "ref_len", "qry_len", 
                  "identiy", "ref_tag","qry_tag" )
```

获取x和y轴的范围

```r
x_range <- range(c(df$ref_start, df$ref_end))
y_range <- range(c(df$qry_start, df$qry_end))
```

新建一个画图设备，并设置好画图区域

```r

plot.new()

plot.window(xlim = x_range,
            ylim = y_range)
```

之后逐行绘制每个联配结果，用不同的颜色来展示倒置的情况

```r
for( i in 1:nrow(df)){
  
  if (df[i,3] < df[i,4]){
    lines(x = df[i,1:2], y = df[i,3:4], col = "red")
  } else{
    lines(x = df[i,1:2], y = df[i,3:4], col = "blue")
  }
  
}
```

最后加上x轴和y轴的标签

```r
box()
axis(1, at = seq(0, x_range[2], 10000), labels = seq(0, x_range[2], 10000) / 10000)
axis(2, at = seq(0, y_range[2], 10000), labels = seq(0, y_range[2], 10000) / 10000)
```

结果如下

![最终结果](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/10/image-eb79ee30ae09459cbe519524a1d63c46.png)

上面这种作图方式就是R语言的基础作图模式，从绘图思想上叫做画笔模式，这区别于ggplot2所代表的图形语法。