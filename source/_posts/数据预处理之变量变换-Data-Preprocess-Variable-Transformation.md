---
title: 数据预处理之变量变换
date: 2019-11-02 18:37:55.904
updated: 2019-11-04 20:22:01.716
url: /archives/Data-Preprocess-Variable-Transformation
categories: 数据科学
tags: 数据挖掘
---

在学习「数据挖掘导论」的**数据预处理**时，里面谈到了变量变换，我联想到了在基因表达量分析时的常见操作，例如FPKM，TPM，CPM，log对数变换。 比如说在文章里面会见到如下的描述

>  The size factor of each cell was computed using a pooling strategy implemented in the R function `computeSumFactors`. Normalized counts were then computed by dividing the counts for each cell by the size factor for that cell. A log2 transformation was applied to normalized counts. 

变量转换有什么好处，需要注意些什么呢？「数据挖掘导论」讨论了两种重要的变量变换类型: 简单函数变换和规范化。

**简单变换**是使用一个简单数学函数分别作用于每一个值，例如log转换，求绝对值，求倒数等。统计学中，变量变换（例如log转换）常用于将不具有高斯（正态）分布的数据变换成具有高斯（正态）分布的数据。在数据挖掘领域可以用来进行数据压缩。这类变换需要我们了解数据在变化前后的后果，例如负数取倒数之后的大小关系会逆转。

**标准化(standardization)或规范化(normalization)**的**目标**是使整个值的集合具有**特定的属性**。使用这两个术语需要特别注意使用这两个词的上下文。书中提到"在数据挖掘界，这两个术语常常可互换，然而，在统计学中，术语规范化可能与使得变量**正态**化相混淆"。虽然从中文的角度来看，规范化和正态化明显不一样，但是从英文的角度看，正态分布翻译自normal distribution, 很容易从normal联想到normalization。当然也有将normalization翻译成归一化，即将原来数据值通过$(x - min(X) / (max(X) - min(X))$转换到0,1之间。统计学的变量标准化指的就是对原来的数据基于均值和标准差计算z-score(公式如下), 不过考虑到均值和标准差受到离群点波动很大，可以用中位数替代均值，用绝对标准差替代标准差。

$$
x' = \frac{x-\bar{x}}{\sigma}
$$


R语言中做标准化常用到一个函数`scale`，它的功能是对矩阵的列进行中心化和(或)缩放

```bash
scale(x, center = TRUE, scale = TRUE)
```

参数就3个，基本上我用的时候就提供第一个输入，举个例子

```r
> x <- c(1,2,3,4,5)
> y <- scale(x)
> y
           [,1]
[1,] -1.2649111
[2,] -0.6324555
[3,]  0.0000000
[4,]  0.6324555
[5,]  1.2649111
attr(,"scaled:center")
[1] 3
attr(,"scaled:scale")
[1] 1.581139
```

之前我都是这样子用的，也没有思考过这到底它到底做了些什么。最近为了彻底搞懂它的两个参数，我翻了它的源代码。

代码总共38行，第一部分是关于`center`参数，也就是中心化。它先判断`center`是否是逻辑值还是数值向量。如果是逻辑值，则当`center=TRUE`时它会计算每一列的均值，然后调用`sweep`按列分别减去均值。如果判定`center`是一个数值向量时，就是为`sweep`手动指定每一列的中心值。

```r
if (is.logical(center)) {
  if (center) {
    center <- colMeans(x, na.rm=TRUE)
    x <- sweep(x, 2L, center, check.margin=FALSE)
  }
}
else {
  if(!is.numeric(center)) center <- as.numeric(center)
  if (length(center) == nc)
    x <- sweep(x, 2L, center, check.margin=FALSE)
  else
    stop("length of 'center' must equal the number of columns of 'x'")
}
```

第二部分是关于`scale`,决定数据如何缩放。如果是逻辑值且为真时，默认用下面的函数计算每一列的缩放系数，然后用`sweep`按列除以每列系数

$$
\sqrt{\frac{\sum_{i=1}^{n}x{'}_i^2}{n-1}}
$$

> 这个计算的结果就是样本标准差，因为$x'_i=x_i-\bar x$

如果是数值向量，则用`sweep`以给定数值按列相除。

```r
if (is.logical(scale)) {
  if (scale) {
    f <- function(v) {
      v <- v[!is.na(v)]
      sqrt(sum(v^2) / max(1, length(v) - 1L))
    }
    scale <- apply(x, 2L, f)
    x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
  }
}
else {
  if(!is.numeric(scale)) scale <- as.numeric(scale)
  if (length(scale) == nc)
    x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
  else
    stop("length of 'scale' must equal the number of columns of 'x'")
}
```

最后一部分是设置输出结果的属性。

```r
if(is.numeric(center)) attr(x, "scaled:center") <- center
if(is.numeric(scale )) attr(x, "scaled:scale" ) <- scale
```

那么默认参数下，我们就是对矩阵按列进行z-score的标准化。检验标准很简单，计算标准化的数据的均值和标准差，因为z-score的结果就是通过线性变换让数据的均值为0，标准差为1（不会改变原来的数据分布）。

```r
> mean(y)
0
> sd(y)
1
```

这解决了我多年的疑惑，为啥我在R语言中就是找不到`zscore`这个函数，因为`scale`默认情况下就是实现`z-score`的变换。当然R语言取名`scale`也是有它的道理，比方说数据中存在离群值，我们应该选择不容易受离群值影响的中位数替代均值，选择绝对平均偏差, 中位数绝对偏差或四分位数极差来替换方差。

```R
# 默认参数
> x <- matrix(c(1,2,3,4,5,1,2,3,4,1000), ncol=2)
> scale(x)

           [,1]       [,2]
[1,] -1.2649111 -0.4505747
[2,] -0.6324555 -0.4483330
[3,]  0.0000000 -0.4460914
[4,]  0.6324555 -0.4438497
[5,]  1.2649111  1.7888488
...
# 替换后
> x <- matrix(c(1,2,3,4,5,1,2,3,4,1000), ncol=2)
> scale(x, 
      center = apply(x, 2, median),
      scale = apply(x, 2, function(x){
        sum(abs(x - mean(x))) / (length(x)-1)
      })
      )

           [,1]         [,2]
[1,] -1.3333333 -0.005012531
[2,] -0.6666667 -0.002506266
[3,]  0.0000000  0.000000000
[4,]  0.6666667  0.002506266
[5,]  1.3333333  2.498746867
...
```

此外，它还能实现其他形式的数据转换，例如归一化，即将数据的范围缩放到0到1之间

```r
x <- matrix(c(1,2,3,4,5,6), ncol=2)
y <- scale(x, center = apply(x, 2, min), 
           scale = apply(x,2,function(x) {max(x)-min(x)})
           )
```

最后总结一下:

- 数据的变量转换有两类，简单变换和标准化/规范化，无论是何种变换都要看它的具体公式
- R语言的`scale`默认就是计算原始数据的z-score, 通过调整`center`和`scale`可以实现多种形式的数据转换。

