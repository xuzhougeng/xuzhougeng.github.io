---
title: 给R使用者的C++最少必要知识
date: 2019-10-13 14:49:10.474
updated: 2019-10-16 09:37:47.716
url: /archives/C++_For_R_User
categories: R
tags: C/C++
---


`C++`是一门非常复杂的编程语言，但是如果你已经有一定的R语言基础，希望通过`C++`来提高代码效率，那么掌握下面几点就能够开始写代码了，然后通过Rcpp调用。

## `C++`是静态编程语言

`C++`是一门静态编程语言，这意味着对于一个变量而言，你需要**先声明**它，才能调用它。而且这个变量名的数据类型在使用过程中是无法更改的，这是因为它在内存中的大小已经固定了。 

举个例子，在R语言中，下面这个代码运行时不会报错（但是会被吐槽）

```r
a <- 1
a <- 'a'
```

但是在`C++`中，下面的代码会被提示`error: redefinition of 'x' with a different type: 'char' vs 'int'`

```C++
char a = 'a';
int a = 0;
```

**注**：每个`C++`语句后需要有`;`, 而R没有。

关于数据结构，标准`C++`的数据结构有如下几种

- 整型: int, long
- 浮点型: float, double
- 逻辑: bool
- 字符: char

`C++`的字符串在STL(Standard Template Library)中。STL里有很多高级数据结构，例如向量(vector)和列表(List)限于篇幅请自行检索。

## `C++`的控制结构

`C++`的控制结构包括，for循环，while循环，if条件语句，switch条件语句，代码形式如下

```C++
// for loop
for ( int i = 1; i <= 10; i++){
    ...;  
}

// while loop
while ( i <= 10){
  ...;
}

// if
if (){
  ...;
} else if(){
  ...;
} else{
  ...;
}

// switch
char grade = 'D';
switch (grade ){
  case 'A':
    cout << 'great' << endl;
    break ;
  case 'B' :
    cout << 'not bad' << endl;
    break ;
  case 'D' :
    cout << 'Not good' << endl;
  default:
    cout << 'wrong' << endl;
    
}
```

在代码形式上，if和while都是和R语言类似，for和switch则是有些不同。

`C++`的循环效率远远高于R语言，而且将对应的R代码修改成`C++`并不复杂。因此会用`C++`的循环，就能优化很多R代码。

## 函数

`C++`的函数和R不同，它的写法如下

```C++
数据结构 函数名(数据结构 参数, ...){
  代码;
  return 变量名;
}
```

举个例子

```C++
int sumTwo(int x, int y){
  int z = x + y;
  return z;
}
```

如果是在R中则是

```R
sumTwo <- function(x, y){
   z <- x + y
   return(z)
}
```

因此在改成原有的R函数时要注意两者的区别，

- 代码形式不同
- 参数名要声明数据类型
- return后直接跟变量名

在函数的使用上，`C++`和R就没有区别了，都是`函数名(变量)`，举个例子`sumTwo(1,2)`

## 额外库加载

在使用R语言的大部分时间里，我们都是面向R包编程，也就是搜索一个R包，然后安装R包，加载R包，调用函数。对于`C++`而言，也有许多现成的库，能够让我们避免造轮子。

举个例子，加载Rcpp库，调用R函数

```C++
#include <Rcpp.h>
#include <cstdio>

using namespace Rcpp;

//[[Rcpp::export]]
int main(){
  printf("N(0,1) 95th percential %9.8f\n",
         R::qnorm(0.95, 0.0, 1.0, 1, 0));

}
```

这里`#include`类似于安装R包，而`using namespace Rcpp;`则是`library`加载R包。

## 面向对象

`C++`和R语言都是一门面向对象编程语言。R里面有S3，S4，RC等形式，`C++`则是`struct`和`class`两种形式，举个例子

```C++
struct Date {
  unsigned int year;
  unsigned int month;
  unsigned int date;
};

class Date {
private:
  unsigned int year;
  unsigned int month;
  unsigned int date;
public:
  void setDate(int y, int m, int d);
  int getDay();
  int getMonth();
  int getYear();
}
```

强行联系的话，R的S3和`C++`的struct比较像，R的S4和`C++`的Date比较像，因为R的S3比较宽松，是基于list的堆叠，内部数据直接暴露到外部，而S4则是比较安全。`C++`的class用private保证private无法直接被外部访问，只能功过public暴露的函数进行操作。

## 指针和内存管理

指针和内存管理是两个比较高级的话题，如果你学习C语言不懂指针那你就和没学过C一样。不过我们学的是`C++`，只要不涉及到很高级的操作，那么我们完全可以避免接触它们。

## 总结

对于一个R语言使用者而言，如果想学习`C++`，那么至少需要知道如下几个概念

- 变量声明
- 控制结构，for, while, if-else
- 函数编写和调用
- 额外库
- 面向对象










