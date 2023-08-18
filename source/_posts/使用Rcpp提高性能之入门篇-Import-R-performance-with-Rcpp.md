---
title: 使用Rcpp提高性能之入门篇
date: 2019-10-15 22:11:19.46
updated: 2019-10-15 22:11:19.46
url: /archives/Import-R-performance-with-Rcpp
categories: R
tags: C/C++
---

`C++`能解决的瓶颈问题有:

- 由于迭代依赖于之前结果，循环难以简便的向量化运算
- 递归函数，或者是需要对同一个函数运算成千上万次
- R语言缺少一些高级数据结构和算法

我们只需要在代码中写一部分`C++`代码来就可以处理上面这些问题。后续操作在Windows下进行，你需要安装Rtools，用`install.packages("Rcpp")`安装新版的Rcpp，最重要一点，你需要保证你R语言时不能是`C:/Program Files/R/R-3.5.1/`这种形式，否则会报错。

后续操作会用到`microbenchmark`包来评估R代码和RCPP的效率差异，用`install.packages('microbenchmark)`安装

## RCPP入门

先从一个简单的`add`函数开始，学习如何用`cppFunction`在R里面写`C++`代码

```r
library(Rcpp)

cppFunction('int add(int x, int y, int z) {
  int sum = x + y + z;
  return sum;
}')
add
# function (x, y) 
# .Call(<pointer: 0x0000000063c015a0>, x, y)
```

Rcpp将会编译`C++`代码， 然后构建能够连接到C++函数的R函数。后续将会介绍如何将一些R代码改写成`C++`代码。

- 标量输入，标量输出
- 向量输入，标量输出
- 向量输入，向量输出
- 矩阵输入，向量输出

### 没有输入，标量输出

最简单的函数就是不提供任何输出，返回一个输出，比如说

```r
one <- function() 1L
```

等价的`C++`代码是

```cpp
int one(){
    return 1;
}
```

那么将这段`C++`代码在R用`cppFunction`中改写就是如下

```r
cppFunction('int one(){
  return 1;
}')
```

上面这段函数就展示了R和`C++`之间一些重要区别:

- `C++`写代码不是` 函数名 <- function(参数){}` 而是 `函数名(函数参数){}`
- `C++`中必须声明返回类型，`ini`就是标量整数。C++对应R语言常用向量的类是: `NumericVector`,`IntegerVector`, `CharacterVector`和`LogicalVector`.
- R语言没有标量，全是向量。而`C++`有向量和标量之分，标量的数据类型是`double`, `int`, `String`和`bool`
- `C++`你必须要用到`return`声明要返回的数据
- 每段代码后要跟着`;`

### 标量输入，标量输出

我们可以写一个函数，`sign`，他的功能就是把一个负数转成正数，正数不变

```r
signR <- function(x){
  if (x > 0){
    x
  } else if (x == 0 ){
    0
  } else{
    -x
  }
}

cppFunction('int signC(int x){
            if( x >0 ){
              return x;
            } else if (x == 0){
              return 0;
            } else {
              return -x;
            }
}')
```

这个例子中要注意两件事情

- 在`C++`中，你需要声明输入的数据类型
- `C++`和R的条件语句长得一样。

### 向量输入，标量输出

R和`C++`一大区别就是R的循环效率很低。因此在R语言要尽量避免使用显示的循环语句，尽量向量化运算函数。而`C++`的循环花销特别小，所以可以放心大胆的用。

让我们用R代码写一个求和函数`sum` 以及 `C++`的求和函数，然后比较下效率

```r
sumR <- function(x){
  total <- 0
  for (i in seq_along(x)){
    total <- total + x[i]
  }
  total
}

cppFunction('int sumC(NumericVector x ){ 
            int n = x.size();
            double total = 0;
            for(int i = 0; i < n; ++i){
              total += x[i];  
            }
            return total;
            }')
```

`C++`版本和R版本的逻辑相同，但是有如下不同

- 用`.size()`确认向量的长度
- `for`的写法为`for(初始值; 判断语句; 递增)`
- **记住**: `C++`的向量索引从0开始，R是从1开始
- 向量赋值是`=`而不是`<-`
- `total += x[i]`等价于`total = total + x[i]`, 类似的符号还有`-=`, `*=`, `/=`

最后用`microbenchmark`比较下,R自带求和函数和我们自己写的两个版本的差异

```r
x <- runif(1000)
microbenchmark(
  sum(x),
  sumC(x),
  sumR(x)
)
```

最快的是高度优化过的内置函数，最差的就是`sumR()`, 速度会比`sumC()`慢10倍以上。 

### 向量输入，向量输出

R中比较常见的操作就是向量间运算，尤其R还会自动补齐。自动补齐某些时候会造成一些问题，但是C++不存在这个问题。我们可以写一个RCPP的`+`函数

```r
cppFunction('NumericVector addC(NumericVector x, NumericVector y){
  int xn = x.size();
  int yn = y.size();
  
  if (xn != yn){
    stop("input should be same length");
  }
  NumericVector out(xn);
  for(int i=0; i< xn; ++i){
    out[i] = x[i] + y[i];
  }
  return out;
}')

x <- runif(1e6)
y <- runif(1e6)
microbenchmark(addC(x,y),
               x+y)
```

### 矩阵输入，向量输出

每个向量类型都有矩阵等价类，`NumericMatrix`, `IntegerMatirx`, `CharacterMatirx`, `LogicalMatirx`. 让我们尝试写一个`rowSums()`函数

```r
cppFunction('NumericVector rowSumsC(NumericMatrix x){
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);
  
  for(int i = 0; i < nrow; i++){
    double total =0;
    for(int j =0; j< ncol; j++){
      total += x(i,j);
    }
    out[i] = total;
  }
  return out;
}')
set.seed(1024)
x <- matrix(sample(100), nrow = 10)
rowSumsC(x)
```

这里注意有两点不同，在`C++`中，你用`()`对矩阵取值，而不是`[]`


尽管看起来`C++`的代码运行起来比R语言快多了，比如说R要一分钟，RCPP只要一秒，但是如果算上我们写代码的时间和调试代码的时间，刚开始不熟练估计要10分钟，那么总体来看，还是直接上手写R代码比较合适。

但是如果有一些代码要不断复用，那么写`C++`代码还是很划算。这个时候就建议将代码写到专门的文本中，用`sourceCpp()`加载，而不是`cppFunction()`函数

在Rsutdio中可以创建一个`C++`模板文件，代码写完之后还可以进行debug。

![创建模板](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/10/image-17ff0fcaeda34cf5bbc0c7676ee3aa11.png)

比如说在里面写上面的`rowSumsC`函数，分为如下几个部分

导入头文件，加载Rcpp到命名空间中，类似于`library()`

```cpp
#include <Rcpp.h>
using namespace Rcpp;
```

使用`// [[Rcpp::export]]`说明这里的函数会被R使用

```cpp
// [[Rcpp::export]]
NumericVector rowSumsC(NumericMatrix x){
  int ncol = x.ncol(), nrow = x.nrow();
  NumericVector out(nrow);
  
  for (int i =0; i < nrow; i++ ){
    double total = 0;
    for (int j =0 ;j < ncol; j++){
      total += x(i,j);
    }
    out[i] = total;
  }
  return out;
}
```

下面部分会在`sourceCpp()`加载后自动运行

```cpp
/*** R
library(microbenchmark)
set.seed(1014)
x <- matrix(sample(100), 10)
microbenchmark(
  rowSumsC(x),
  Matrix::rowSums(x)
)
*/
```

将文件保存成`rowSumsC.cpp`, 之后在R里用`sourceCpp(file = "rowSumsC.cpp")`。

