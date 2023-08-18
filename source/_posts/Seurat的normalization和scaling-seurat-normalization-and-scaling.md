---
title: Seurat的normalization和scaling
date: 2019-12-17 16:07:37.785
updated: 2019-12-17 16:07:43.354
url: /archives/seurat-normalization-and-scaling
categories: R
tags: Seurat | 单细胞 | 数据挖掘
---


Seurat的分析流程有两步, 对数据的normalization和scaling. 两种的作用不同，前者是为了处理每个细胞的总count不同的问题，而后者则是让每个基因的表达量的均值为0，方差为1.

normlization对应的函数是`NormalizeData`，通过数据进行一些列变换，消除**文库大小**的影响。 它有三种方法, LogNormalize, CLR, RC

默认方法方式是**LogNormalize**， 即对于每个细胞，将每个基因的count除以总数，然后乘以一个scale.factor, 之后以自然对数进行转换。为了提高效率，Seurat编写了`C++`代码用于加速

```c++
// [[Rcpp::export]]
Eigen::SparseMatrix<double> LogNorm(Eigen::SparseMatrix<double> data, int scale_factor, bool display_progress = true){
  Progress p(data.outerSize(), display_progress);
  Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      it.valueRef() = log1p(double(it.value()) / colSums[k] * scale_factor);
    }
  }
  return data;
}
```

也就是说，如果你不追求效率，我们是可以用纯R代码实现。

```r
mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
mat_norm <- LogNormalize(data = mat)
# LogNormalize等价于
log1p(t(t(mat) / colSums(mat)) * 10000)
```

Scale这一步对应的函数是`ScaleData`， 在它处理之后，使得每个基因在所有样本的均值是0，而方差是1。这让每个基因在下游分析中的具有相同的权重，使得高表达基因不那么显著。

而如果需要移除数据集中不需要的变异来源(unwanted sources of variation),  `ScaleData`需要设置额外的参数`vars.to.regress`，但是作者更加推荐使用`SCTransform`。

尽管`ScaleData`对应的源代码非常的长，但是和数据处理相关的就是如下这几条

```r
# Currently, RegressOutMatrix will do nothing if latent.data = NULL
data.scale <- scale.function(
  mat = object[features[block[1]:block[2]], split.cells[[group]], drop = FALSE],
  scale = do.scale,
  center = do.center,
  scale_max = scale.max,
  display_progress = FALSE
)
```

其中`scale.function`是根据数据类型来决定是`FastRowScale`还是`FastSparseRowScale`, 而这两个代码Seurat也是通过`C++`进行提速了。因此，如果不追求效率，还是可以用纯R代码实现。

```r
pbmc_small <- ScaleData(pbmc_small)
GetAssayData(pbmc_small, "scale.data")[c("IGLL5"),1:2]
# 等价于
scaled_data <- t(scale(t(GetAssayData(pbmc_small, "data"))))
scaled_data[c("IGLL5"),1:2]
```

> 默认R语言的`scale`是按列处理，对于行为基因，列为样本的数据，可以直接套用。但是Seurat的输入数据是行为基因，列为样本，那么就需要按行scale。换句话说，如果你需要对行进行scale，那么你可以通过`Seurat:::FastRowScale()`的方式调用Seurat写的`FastRowScale`函数。此外Seurat编写了大量的`Fast*`函数，都可以尝试用在代码中。


