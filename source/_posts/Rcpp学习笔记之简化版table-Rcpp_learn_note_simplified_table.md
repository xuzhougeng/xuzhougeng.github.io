---
title: Rcpp学习笔记之简化版table
date: 2019-10-14 13:53:46.564
updated: 2019-10-14 13:57:30.458
url: /archives/Rcpp_learn_note_simplified_table
categories: R
tags: C/C++
---


R语言有一个自带的函数`table`能够统计输入变量中不同元素出现的次数，举个例子

```r
d <- rep(c("A","B","C"), 10)
table(d)
```

果子老师曾写一篇推送，自己写了一个简化版的table，比R自带的table 运行的速度更快，如下

```r
tableGZ <- function(x){
  if(sum(is.na(x)) == 0){
    data <- x
    input <- unique(x, fromLast = TRUE)
    dd <- sapply(input, 
                 function(x) {sum(data==x)})
    names(dd) <- unique(data, fromLast = TRUE)
    dd
  } else{
    data <- x[!is.na(x)]
    input <- unique(x, fromLast = TRUE)
    dd <- sapply(input, function(x){
      sum(data == x)
    })
    dd <- c(dd, sum(is.na(x)))
    names(dd) <- c(input, 'NA')
    dd
  }
}
```

我们通过运行1000次代码，来比较下两者的运行速度

```r
bench::system_time(for ( i in 1:1000){
  tableGZ(d)
})
# 结果
process    real 
 42.9ms  42.4ms 

bench::system_time(for ( i in 1:1000){
  table(d)
})
# 结果
process    real 
  107ms   106ms  
```

在我的电脑上，果子老师的代码运行速度比R自带的table快了将近3倍。当然这是有原因的，因为R的table的代码功能更加复杂，能够比较多个变量之间的关系，例如`table(d,d)`。

既然是简单的统计每个元素的次数，那么我就想着能不能写出一个比果子老师速度更快的函数。 于是，我抽空写了下面的代码


```r
tableZG <- function(x){
  
  NA_pos <- is.na(x)
  NA_num <- sum(NA_pos)
  
  x <- x[!NA_pos]
  
  out <- vector(length = length(x))
  out_name <- rep(NA,  length(x))
  
  for (i in 1:length(x)) {
    
    for (j in 1:length(out_name)){
      if ( is.na(out_name[j]) ){
          out_name[j] <- x[i]
          out[j] <- 1
          break
        }
      
      if ( out_name[j] == x[i] ){
        out[j] <- out[j] + 1
        break
      } 
    }
  }
  
  if (NA_num > 0){
    na_end <- sum(!is.na(out_name))
    out_name[na_end + 1] <- 'NA'
    out[na_end + 1] <- NA_num
    
  } 
  na_pos <- is.na(out_name)
  out_name <- out_name[!na_pos]
  out <- out[!na_pos]
  names(out) <- out_name

  return(out)
  
}
```

虽然我的代码更长了，但是并没有让速度提高，反而比果子老师的代码慢，甚至还不如R自带的table。

```r
system.time(for ( i in 1:1000){
  tableZG(d)
})
# 结果
process    real 
  148ms   148ms 
```

当然那么一长串代码并不是白写的，因为我特意避免了使用R特有的内容，所以代码能够很容易改写成`C++`代码使用`Rcpp`调用，从而提高速度

新建一个`tableC.cpp`文件，代码内容如下

```C++
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector tableC(CharacterVector cv){
  
  // initialize variable
  CharacterVector na = CharacterVector::create("NA");
  NumericVector out = rep(NumericVector::create(0), cv.size());
  CharacterVector out_name = rep(na, cv.size());
  int unique_num = 0;
  
  // 
  for (int i = 0; i < cv.size();i ++) {
    
    for (int j = 0; j < cv.size(); j++){
      
      if ( out_name[j] == "NA" ){
        out_name[j] = cv[i] ;
        out[j] = 1;
        unique_num += 1;
        break;
      }
      
      if ( out_name[j] == cv[i] ){
        out[j] = out[j] + 1;
        break;
      } 
    }
  }
 
  out.attr("names") = out_name;
  
  return out;
  
} 
```

然后在R里面用Rcpp这个C++代码，替换掉之前代码中的循环部分

```r
Rcpp::sourceCpp("tableC.cpp")
tableZG <- function(x){
  
  NA_pos <- is.na(x)
  NA_num <- sum(NA_pos)
  
  x <- as.character(x[!NA_pos])
  res <- tableC(x)
  res <- res[!names(res) == "NA"]
  
  if (NA_num > 0){
    res <- c(res, "NA"=NA_num)
  }
  return(res)
}
```

于是这一次在C++的加持下，我写的table函数速度超过了果子老师的代码。

```r
bench::system_time(for ( i in 1:1000){
  tableZG(d)
})
#结果
process    real 
 30.1ms  29.3ms 
```

最后总结一下：如果一个操作只需要做一次，那么速度可能并不是最重要的。因为即便是一个原本要花24小时的代码，提速10倍，只要2小时，你可能也会愿意等一等。但是如果这个操作需要重复很多次，上百次，上千次，甚至上万次，那么你就可能等不下去了。你就需要对代码中的一些限速步骤进行优化，比如说table这种多功能函数，你就可以自己用R写一个简化版的函数，替换掉原先的代码。

如果对速度有更高的要求，那么就需要用到`C++`进行代码重写了。学习`C++`其实并不会特别难，因为有一个`Rcpp`简化了许多操作，你只需要掌握几个最基本的语言特性，比如说`C++`需要先定义变量才能使用变量。

## 推荐阅读

- [给R使用者的C++最少必要知识](/archives/C++_For_R_User)
- [Rcpp学习笔记之Hello World!](/archives/Write_first_function_Using_Rcpp)