---
title: 使用mclust进行聚类分析
date: 2019-11-09 16:05:10.766
updated: 2019-11-09 16:24:59.416
url: /archives/Cluster-data-with-mclust
categories: R
tags: 数据挖掘
---


`mclust`(Model-based clustering) 能够基于高斯有限混合模型进行聚类，分类以及密度估计(density estimation)。对于具有各种协方差结构的高斯混合模型，它提供了根据EM算法的参数预测函数。它也提供了根据模型进行模拟的函数。还提供了一类函数，整合了基于模型的层次聚类，混合估计的EM算法，用于聚类、密度估计和判别分析中综合性策略的贝叶斯信息判别标准。最后还有一类函数能够对聚类，分类和密度估计结果中的拟合模型进行可视化展示。

简而言之，`mclust`在R语言上实现了**基于高斯有限混合模型**的**聚类**，**分类**和**密度估计**分析，并且还有专门的可视化函数展示分析结果。

> 和`mclust`功能相似的其他R包: 'Rmixmod', 'mixture', 'EMCluster', 'mixtools', 'bgmm', 'flexmix'

## 安装和加载

在已有的R语言的基础上，只需要运行如下代码即可

```R
# 安装
install.packages("mclust")
# 加载
library(mclust)
```

## 聚类实战

以一个例子来介绍一下如何使用`mclust`进行聚类分析。我们得要先加载一个来自于R包'gclus'的数据集'wine'，该数据集有178行，分别是不同区域的品种葡萄， 14列，其中后13列是化学分析的测量值。我们的目标是将其进行分类。

第一步: 加载数据集

```r
install.packages("gclus")
data("wine", package = "gclus")
dim(wine)
```

第二步 : 使用`Mclust`做聚类分析. `Mclust`主要功能就是分析当前的提供的数据是由什么统计模型

```r
# 第一列和聚类无关
X <- data.matrix(wine[,-1])
mod <- Mclust(X)
```

直接在交互行输入`mod`会得到如下信息

```r
'Mclust' model object: (VVE,3) 

Available components: 
 [1] "call"           "data"           "modelName"     
 [4] "n"              "d"              "G"             
 [7] "BIC"            "bic"            "loglik"        
[10] "df"             "hypvol"         "parameters"    
[13] "z"              "classification" "uncertainty" 
```

这里需要对结果稍作解释，第一行告诉我们'Mclust'以`VVE`模型将数据分为3类。第3行开始，它告诉我们'Mclust'的输出结果中包含了如下内容，我们可以通过`$`来提取。举个例子，我们提取`Mclust`的聚类结果和已知结果进行比较

```r
table(wine$Class, mod$classification)
# 如下是输出信息
     1  2  3
  1 59  0  0
  2  0 69  2
  3  0  0 48
# adjustedRandIndex:评估聚类效果
adjustedRandIndex(wine$Class, mod$classification)
```

从结果中，我们发现仅有2例没有正确聚类，说明`Mclust`的效果很好。但是随之而来的问题是，`Mclust`如何挑选模型以及它为什么认为聚成3类比较合适呢？我们可以根据什么信息进行模型选择呢？

## 模型选择

为了解答上面的问题，我们需要稍微了解点`Mclust`的原理。和其他基于模型的方法类似，`Mclust`假设观测数据是一个或多个混合高斯分布的抽样结果，`Mclust`就需要根据现有数据去推断最优可能的模型参数，以及是由 q几组分布抽样而成。`mclust`一共提供了14种模型（见下表），可以用`?mclustModelNames`可以查看`mclust`提供的所有模型。

![模型总结](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/11/image-20191109125815419-9f69302312944a0f92e0c245c6a39f0f.png)

为了对模型有一个直观的理解，`mclust`提供了这些模型数据分为三组前提下在二维中的形状。

![二维的模型形状](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/11/image-20191109125852039-00b71d5c1283487c8084b329c61fb7aa.png)

继续回到之前的问题，`Mclust`如何确定模型和确定分组数目。之前我们调用`Mclust`时，除了必须设置的输入参数，没有修改其他参数。其实`Mclust`可以设置的参数不少，和问题直接相关的是如下两个参数

- G: 分组数，默认情况下是`1:9`
- modelNames: 待拟合的模型，默认使用所有14种。

也就是，`Mclust`默认得到14种模型1到9组的分析结果，然后根据一定的标准选择最终的模型和分组数。

`Mclust`提供了两种方法用于评估不同模型在不同分组下的可能性

- BIC( Bayesian Information Criterion ): 贝叶斯信息判别标准
- ICL( integrated complete-data likelihood ): 综合完全数据可能性

`Mclust`默认用的就是`BIC`，因此我们可以用`plot.Mclust`绘制其中BIC变化曲线

```r
plot.Mclust(mod, what = "BIC", 
     ylim = range(mod$BIC[,-(1:2)], na.rm = TRUE), 
     legendArgs = list(x = "bottomleft", cex =0.7))
```

![BIC曲线](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/11/image-20191109154734889-6583692493584d749ae401eab1d6889a.png)

`Mclucst`会选择其中BIC最大的模型和分组作为最终的结果。

此外我们可以用`MclustBIC`和`MclustICL`分别进行计算

```r
par(mfrow=c(1,2))
BIC <- mclustBIC(X)
ICL <- mclustICL(X)
```

从中选择最佳的模型分组和模型作为输入

```r
mod2 <- Mclust(X, G = 3, modelNames = "VVE", x=BIC)
```

## 可视化展示

`mclust`为不同的输出都提供了对应的泛型函数用于可视化，你需要用`plot`就能得到结果。例如对之前的聚类结果在二维空间展示

```r
drmod <- MclustDR(mod, lambda = 1)
plot(drmod)
# 会提供一些列选项让你选择, 展示不同的结果
# Dimension reduction for model-based clustering and classification plots: 

1: scatterplot
2: contour
3: classification
4: boundaries
5: density
6: evalues
# 以1为例
```

![plot](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/11/image-20191109155624344-40f3e047fc6644c586e932b9196fcce3.png)

`mclust`还有很多其他功能，例如密度估计，自举推断(Bootstrap inference)，这些内容建议阅读"mclust 5: Clustering, Classification and Density Estimation Using Gaussian Finite Mixture Models "


## 推荐阅读

想要更好的学习这个R包的使用，还需要去学习如下概念

- EM算法( expectation–maximization algorithm )
- BIC
- MLE(maximum likelihood estimator)