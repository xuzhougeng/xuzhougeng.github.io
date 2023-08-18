---
title: log化的TPM能做差异分析吗
date: 2022-04-22 02:43:10.256
updated: 2022-04-22 02:43:10.256
url: /archives/differential-analysis-in-log-tpm-data
categories: 数据科学
tags: 转录组
---

有人提了一个问题: 我用log2TPM做差异，log2本来就是对TPM标准化的嘛，我看样品表达分布还是不是特别齐，还需要再进一步标准化吗

对于这个问题，我的回答如下: 

首先，我们如果要对转录组做差异分析，那么最好的数据应该就是raw count，而不是这些转换后的数据。这是因为差异分析的本质是做统计检验，而我们应该都知道统计检验都会对数据有一些基础假设，例如，我们平常用的t-test，默认假设是数据服从正态分布。常用的DESeq2, limma(voom), edgeR之所以都要求输入数据是count，就是因为它们的统计假设是基于count。

其次，如果我们只能拿到TPM或log2TPM的数据，那在这个基础上做差异分析所能用的统计检验方法就建议是无参方法，例如wilcoxon秩和检验。尽管log2转换会让偏倚的数据更符合钟形曲线（正态分布），但你不能在基础上用t-test。 举个例子, 

```R
x <- c(100, 200, 300) 
y <- c(1000, 2000, 3000)

t.test(x,y) # p=0.08787
t.test(log(x), log(y)) # p=0.0071

wilcox.test(x,y) # p=0.1
wilcox.test(log(x), log(y)) # p=0.1

```

在我构建的两个数据中，经log转换，t-test从p>0.05变成了p<0.01, 而秩和检验没有变化。

最后，我们在转录组分析中会用到TPM归一化和log2转换的目的并不是为了差异分析，而是为了后续的可视化。比如说绘制热图，我们就不希望数据里面不要有太高和太低的值，使得整体的颜色极端化。

> In order to test for differential expression, we operate on raw counts and use discrete distributions as described in the previous section on differential expression. However for other downstream analyses – e.g. for visualization or clustering – it might be useful to work with transformed versions of the count data.  
[Analyzing RNA-seq data with DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)


综上，我们用log2TPM做差异应该采用无参的方法，而无参的方法并不考虑样本的分布。


推荐阅读:

- [别再用DEseq2和edgeR进行大样本差异表达基因分析了
](https://kaopubear.top/blog/2022-03-20-donot-use-deseq2-edger-in-human-population-samples/)
