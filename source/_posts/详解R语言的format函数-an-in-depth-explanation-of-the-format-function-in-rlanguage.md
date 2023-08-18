---
title: 详解R语言的format函数
date: 2023-04-01 08:03:26.662
updated: 2023-04-01 08:03:26.662
url: /archives/an-in-depth-explanation-of-the-format-function-in-rlanguage
categories: R
tags: 源码解读
---

先提一个问题，请问下面最后的两个逻辑语句的判断，结果是什么？

```Python

x <- c(1234567.8, 12.12345)
x11 <- format(x, scientific = FALSE)
x12 <- format(x, scientific = TRUE)
x21 <- as.numeric(x11)
x22 <- as.numeric(x12)

x21 == x
x22 == x

```

大家可以稍作思考。

之所以给出开头的问题，还是因为之前讨论的apply函数。当apply应用到data.frame上时，中间会调用泛型函数as.matrix将数据进行转换，如果发现数据中存在非数值型，会调用format对数值型数据进行格式化。

```R
....
  else if (non.numeric) {
    for (j in pseq) {
      # 如果是字符串，不转换
      if (is.character(X[[j]])) 
        next
      # 如果是逻辑值，调用as.character
      else if (is.logical(xj <- X[[j]])) 
        xj <- as.character(xj)
      else {
        # 如果有NA，先确定位置
        miss <- is.na(xj)
        # 如果是因子，调用as.vector，否则调用format
        xj <- if (length(levels(xj))) 
          as.vector(xj)
        else format(xj)
        # 最后把NA放回
        is.na(xj) <- miss
      }
      X[[j]] <- xj
    }
  }
....
```

这个format的描述只有一句，‘Format an **R** object for pretty printing.’即，格式化一个 R 对象用于**漂亮**的打印，

问题出在这个“漂亮”上，R语言有自己的理解。

比如说同样的 1234567.8在不同的环境下，会出现了不同的情况，

```Python
format(c(1234567.8, 12.12345))
# "1.234568e+06" "1.212345e+01"
format(c(1234567.8, 123.12345))
# "1234567.8000" "    123.1235"

```

一个采取了科学计数法，并保留了7位有效数字，显示为1.234568e+06，一个则保持原状1234567.8000。

为什么会出现这个情况呢？这和scientific参数有关，接受一个逻辑值或者一个整数。如果提供的逻辑值，行为非常简单，就是要么用科学计数法，要么不用科学计数法。

而问题出在它默认值**不是逻辑值**，而是**NA**，此时他使用 `getOption("scipen")` 的结果作为输入，用于决定是以固定，还是以指数符号（科学计数法）打印数字。正值偏向于固定符号，负值则偏向于科学计数法。

具体这个参数会触发R哪一行底代码呢？通过我刨根问底，最终确定源代码中src/mian/format.c 的formatReal函数，一个用于对double类型数据进行格式化的函数。

代码有一段注释，含义是，在满足有效位数的前提下，**固定表示**使用的空间不超过**指数表示**时，优先使用固定显示。

```C
    /* F Format: use "F" format WHENEVER we use not more space than 'E'
     *    and still satisfy 'R_print.digits' {but as if nsmall==0 !}
     *
     * E Format has the form   [S]X[.XXX]E+XX[X]
     *
     * This is indicated by setting *e to non-zero (usually 1)
     * If the additional exponent digit is required *e is set to 2
     */
     ...
  *d = mxns - 1;
  *w = neg + (*d > 0) + *d + 4 + *e; /* width for E format */
  if (wF <= *w + R_print.scipen) { /* Fixpoint if it needs less space */
      *e = 0;
      if (nsmall > rgt) {
    rgt = nsmall;
    wF = mxsl + rgt + (rgt != 0);
      }
      *d = rgt;
      *w = wF;
  } /* else : "E" Exponential format -- all done above */
...       
```

代码中w是指数表示法的宽度, wF是固定表示的宽度，最后w会加上 scipen 得到最终的长度，默认是0，也就是公平决斗。

回到数据`c(1234567.8, 12.12345)`中，我们来看它在固定表示和指数表示的结果。

```R
> format(c(1234567.8, 12.12345), scientific = TRUE)
[1] "1.234568e+06"
[2] "1.212345e+01"
> format(c(1234567.8, 12.12345), scientific = FALSE)
[1] "1234567.80000"
[2] "     12.12345"

```

在保证最小数字有7位有效数字的前提下（由digits参数控制），采用指数显示法时结果是12个字符，采取固定显示结果是13个字符。当scipen=0时，指数显示**比较窄**，因此选择指数显示。

而数据`c(1234567.8, 123.12345)`，在保证最小数字有7位有效数字时，结果都是12个字符，势均力敌的情况下，采取固定显示。

```R
> format(c(1234567.8, 123.12345), scientific = TRUE)
[1] "1.234568e+06"
[2] "1.231235e+02"
> format(c(1234567.8, 123.12345),scientific = FALSE)
[1] "1234567.8000"
[2] "    123.1235"
```

> 有效数字的定义是，从第一位非0数字开始到最后一位，一共有多少位数字。例如100是3为有效数字，0.001是1位有效数字，0.00100是3位有效数字。

最后帮助大家梳理下整体逻辑

1. format函数作用是让打印结果更加美观
2. 在打印数字时，在保证最小数以一定有效数值展示的前提下，对比固定表示和指数表示两者的字符数目，选择占用字符最少的方法
3. 最终数值的展示效果会受到digits和scientific的影响
4. digits和scientific默认值来自于getOption("digits")和getOption("scipen")
5. 如果想让结果都以固定显示，可以设置scipen为一个非常大的数。

