---
title: chatGPT会替代我们吗？
date: 2022-12-06 07:57:22.472
updated: 2022-12-06 07:57:22.472
url: /archives/is-chatgpt-going-to-replace-us
categories: 其他
tags: 
---

最新上线chatGPT，非常的热闹，openAI的服务器几度因为请求太多而宕机，紧急增加了服务器。很多视频都说程序员岌岌可危，因为chatGPT真的能帮我们写代码。

我试了试，发现chatGPT确实好用，比如说我让他帮我写一段差异分析的代码(differential analysis in R code)

```
# Load the necessary packages
library(DESeq2)

# Import the RNA-seq data into R
data <- readRDS("rnaseq_data.rds")

# Create a DESeqDataSet object from the data
dds <- DESeqDataSet(data, design = ~ condition)

# Perform the differential expression analysis
dds <- DESeq(dds)

# Extract the results of the analysis
results <- results(dds)

# Filter the results to only include significantly differentially expressed genes
results <- results[which(results$padj < 0.05), ]
```

这是他的实现方式，确实有模有样。但其实，我如果去查DESeq2的相关的教程，也能找到类似的代码呀。

另外，有些时候，它还会给出错误的实现方式，比如说，我让他绘制GWAS的Q-Q图，以下是它给的代码。

```
# Load the necessary packages
library(graphics)
library(GWASdata)

# Load the example GWAS data
data(gwas1)

# Create a Q-Q plot for the GWAS data
qqplot(gwas1$pvalue)
```

其中有一个R包，GWASdata需要额外安装。于是，我问他怎么安装呢？他给出的建议是用 `install.packages("GWASdata")`。然而，事实上，GWASdata是放在Bioconductor上的，所以安装方式是`BiocManager::install("GWASdata")`。当然即便我安装上了，但是这个代码也无法运行，因为gwas1这个数据集不存在与GWASdata中，我目前也不知到和这个gwas1哪里来的。更有趣的是，经过我的检索，实际上，我们应该用`qqPlot`这个来自于qqman包的函数，才能根据p值来画图。

实际上，用AI辅佐写代码的事情，之前用GitHub的Copilot时，已经感受到了它的实用性了。它能够帮忙补全很多代码，然后我们在此基础上做些修改，差不多就能用了。但是，真的让它完成整个项目，我觉得还有一段路要走吧。





