---
title: Rcpp学习笔记之Hello World!
date: 2019-10-13 15:52:36.99
updated: 2019-10-13 15:52:36.99
url: /archives/Write_first_function_Using_Rcpp
categories: R
tags: C/C++
---

> 在使用R语言多年以后，我终于开始去学习Rcpp，利用`C++`来提高运行速度。其实当你能熟练的使用一门语言后，再去学一门新的语言，并没有想象中的那么难，更何况Rcpp把很多脏活累活都给包办了，在里面调用`C++`还是挺方便。

学习Rcpp的最重要一步是，运行一个"hello world!" 。如果能够运行"hello world!"就表明搞定了环境配置，后面就可以愉快的写代码了。

安装Rcpp的方式为，`install.packages("Rcpp")`， 安装过程中可能会出现一些问题，对于不同的操作系统需要做不同的准备工作，

- Windows: 安装Rtools
- Mac: 安装Xcode，需要在App store下载
- Linux: 需要有GCC的编译环境

之后就让我们写人生中第一个`C++`函数， `hello`,

```c++
int hello(){
    std::cout << "Hello, World!";
    return 0;
}
```

那么问题来了，这个函数应该如何才能让R语言调用呢？最简单的方式就是Rcpp的`cppFunction`

```r
library("Rcpp")
cppFunction('
int hello(){
    std::cout << "Hello, World!";
    return 0;
}
')
hello()
```

上面语句中，cppFunction中的`''`中的内容是C++语句。在这条语句中，我们以`int hello(){}`定义了一个hello函数，这个函数不接受参数，返回一个整型。在函数里面一共有两条语句，每条语句都以`;`结尾。

第一条是调用了`C++`的标准库的cout, `std::cout`, 和R中以`包名::函数名`调用函数的方法类似。第二条则是`return 0`，返回结果。

除了利用`cppFunction`外，另外一种更常用的方法就是将代码放在其他文件中，然后用`sourceCpp`的方式读取。我们新建一个`hello.cpp`的文件，里面的内容如下(如果用Rstudio新建C++ 文件，它会提供一个模版用于修改)

```c++
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
int hello(){
    cout << "Hello, World!";
    return 0;
}

/*** R
hello()
*/
```

之后在R里面加载并调用，和之前的结果一样。

```r
library(Rcpp)

sourceCpp("hello.cpp")
hello()
```

解释下`C++`的代码。第一行是用`#include <Rcpp.h>`语句导入了Rcpp的头文件，类似于R里面的`install.package()`用于安装R包，然后用`using namespace xxx;` 的方式加载了两个库，Rcpp和std, 这就类似于R里面的`library()`函数。

`C++`里面用`//`和`/*** 注释语句 */`进行代码注释。前者是注释单行，类似于R里面的`#`注释，后者是注释多行。只不过`//[Rcpp::export]]`这条注释有特殊的含义，在`sourceCpp`读取代码解析的过程中，被这条语句注释的函数能够在R里面调用。换句话说，如果你删了这句话，那么这个`hello`函数在R里面就是无法直接调用的。而`/*** R */`里面可以放R代码，会`c++`代码编译结束后运行，常用于代码测试。

假如我们能够成功运行上面的代码，那么接下来要做的事情就是学习`C++`的基本语法(参考[给R使用者的C++最少必要知识](/archives/C++_For_R_User))，学习Rcpp的数据结构。


