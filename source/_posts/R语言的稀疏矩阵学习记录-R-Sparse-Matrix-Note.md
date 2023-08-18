---
title: R语言的稀疏矩阵学习记录
date: 2019-08-17 10:18:28.869
updated: 2019-09-02 12:35:44.438
url: /archives/R-Sparse-Matrix-Note
categories: R
tags: 数据结构
---

一个很大的矩阵， 320127 行,  8189列，假如用一个全为0的普通矩阵来存储，需要用到9.8Gb

```r
cols <- 8189
rows <- 320127
mat <- matrix(data = 0, nrow=320127, ncol = 8189)
print(object.size(mat), unit="GB")
# 19.5 Gb
mat <- matrix(data = 0L, nrow=320127, ncol = 8189)
print(object.size(mat), unit="GB")
# 9.8 Gb这里的0其实也要区分
```

这里的`0L`表示数据类型是`integer`，默认是`numeric`. 这两者最大的区别在于，当你用`320127L * 8189L`，你会得到一个NA，而`320127 * 8189`不会

如果用稀疏矩阵保存的话

```r
mat <- Matrix(data = 0L, nrow=320127, ncol = 8189, sparse = TRUE)
print(object.size(mat), unit="GB")
#0 Gb
dim(mat)
#[1] 320127   8189
```

虽然行列数一样，但是稀疏矩阵几乎不占用任何内存。而且普通矩阵支持的运算，比如说求行和，求列和，提取元素的操作，在稀疏矩阵矩阵也是可以的，只不过会多花一点点时间而已。同时还有很多R包支持稀疏矩阵，比如说`glmnet`，一个做lasso回归的R包。

虽然看起来稀疏矩阵很美好，但是在R语言中那么大的稀疏矩阵的部分操作会出错

```r
> mat2 <- mat + 1
Error in asMethod(object) : 
  Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 105
```

即便是我想把它用`as.matrix`转回普通矩阵，它也报错了

```r
> mat3 <- Matrix::as.matrix(mat)
Error in asMethod(object) : 
  Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 105
```

既然现成的`as.matrix`无法处理，那怎么办呢？最简单粗暴的方法就是新建一个普通矩阵，然后对稀疏矩阵进行遍历，将稀疏矩阵的值挨个放回到的普通矩阵上。

```r
mat2 <- matrix(data = 0, nrow=320127, ncol = 8189)
for (i in seq_len(nrow(mat))){
    for (j in seq_len(ncol(mat))){
        mat2[i,j] <- mat[i,j]
    }
}
```

那么这大概要多少时间呢？反正我的电脑跑了2个小时也没有跑完，所以你也别测试了。

那有没有办法可以加速呢？加速的方法就是减少for循环的次数，因为我们是一个稀疏矩阵，大部分的空间都是0，我们只需要将不为0的部分赋值给新矩阵即可。

这需要我们去了解下稀疏矩阵的数据结构

```r
> str(mat)
Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  ..@ i       : int(0) 
  ..@ p       : int [1:8190] 0 0 0 0 0 0 0 0 0 0 ...
  ..@ Dim     : int [1:2] 320127 8189
  ..@ Dimnames:List of 2
  .. ..$ : NULL
  .. ..$ : NULL
  ..@ x       : num(0) 
  ..@ factors : list()
```

`@Dim`记录矩阵的维度信息, `@Dimnames`记录行名和列名, `@x`记录不为0的数值。`@i`记录不为0的行索引，和`@x`对应，这里全为0，所以不记录。`@p`比较复杂，并不是简单的记录不为0值的列索引，看文档也不知道是啥，不过通过检索可以找到它和不为0值的列索引的换算关系。

因此代码优化为

```r
row_pos <- mat@i+1
col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
val <- mat@x
    
for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
}
```

可以将其封装为一个函数

```r
as_matrix <- function(mat){

  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
    
  for (i in seq_along(val)){
      tmp[row_pos[i],col_pos[i]] <- val[i]
  }
    
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
```

如果速度还需要提高，那么可能就需要Rcpp上场了. 我参考着<http://adv-r.had.co.nz/Rcpp.html>写了一个简单的代码

```r
Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){

  int k = z.size() ;

  IntegerMatrix  mat(nrows, ncols);

  for (int i = 0; i < k; i++){
      mat(rp[i],cp[i]) = z[i];
  }

  return mat;
}
' )


as_matrix <- function(mat){

  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
  
  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                 nrows =  mat@Dim[1], ncols = mat@Dim[2])

  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
```


如果之前的矩阵有78945836个元素，`system.time`显示只需要40s。

