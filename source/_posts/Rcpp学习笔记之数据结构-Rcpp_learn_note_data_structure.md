---
title: Rcpp学习笔记之数据结构
date: 2019-10-15 22:01:05.852
updated: 2019-10-15 22:01:05.852
url: /archives/Rcpp_learn_note_data_structure
categories: R
tags: C/C++
---

> 在使用R语言多年以后，我终于开始去学习Rcpp，利用C++来提高运行速度。其实当你能熟练的使用一门语言后，再去学一门新的语言，并没有想象中的那么难，更何况Rcpp把很多脏活累活都给包办了，在里面调用C++还是挺方便。

C++是一门静态编译面向对象的编程语言，R是动态解释性面向对象语言，那么有一个不同就在于，你需要先声明一个变量，才能调用该变量。而在声明变量的时候，你就会遇到一个R语言中不怎么思考的问题，我这个变量要存放什么样的数据呢？

## 数据类型

R的基本数据类型有六种 "logical", "integer", "numeric" (等价于"double"), "complex", "character" 和 "raw"，但常用的就四种， "integer", "numeric","logical"和"character"。

对于数值数据而言，`C++`的整型和浮点型可以定义多种精度，而R语言则简化成两种"integer"(32 bit)和"numeric"(53 bits)，并且默认情况下数值都是浮点型。

在数据结构上，R语言提供了向量(vector)，矩阵(matrix)，数组(array)，数据框(data.frame)和列表(list)，唯独没有标量。其中向量可以认为是R语言的基本单位，是其他数据结构的基础。

原本能够让`C++`处理R对象，以及将`C++`处理完结果转成R能识别的对象，你需要了解R的内部结构，但是`Rcpp`通过定义了一系列数据类型（如下），使得我们能够非常容易地在R和`C++`之间交换对象。

```C++
// Rcpp/include/Rcpp/vector/instantiation.h
namespace Rcpp{

    typedef Vector<CPLXSXP> ComplexVector ;
    typedef Vector<INTSXP> IntegerVector ;
    typedef Vector<LGLSXP> LogicalVector ;
    typedef Vector<REALSXP> NumericVector ;
    typedef Vector<REALSXP> DoubleVector ;
    typedef Vector<RAWSXP> RawVector ;

    typedef Vector<STRSXP> CharacterVector ;
    typedef Vector<STRSXP> StringVector ;
    typedef Vector<VECSXP> GenericVector ;
    typedef Vector<VECSXP> List ;
    typedef Vector<EXPRSXP> ExpressionVector ;

    typedef Matrix<CPLXSXP> ComplexMatrix ;
    typedef Matrix<INTSXP> IntegerMatrix ;
    typedef Matrix<LGLSXP> LogicalMatrix ;
    typedef Matrix<REALSXP> NumericMatrix ;
    typedef Matrix<RAWSXP> RawMatrix ;

    typedef Matrix<STRSXP> CharacterMatrix ;
    typedef Matrix<STRSXP> StringMatrix ;
    typedef Matrix<VECSXP> GenericMatrix ;
    typedef Matrix<VECSXP> ListMatrix ;
    typedef Matrix<EXPRSXP> ExpressionMatrix ;

}
```

以R语言的向量为例，Rcpp考虑到R语言所有可能的数据类型，提供了9种Vector。比如说"NumericVector", "IntegerVector", "LogicalVector"," CharacterVector"就是对应着"integer", "numeric","logical"和"character"这四种数据类型。因此，你可以根据你输入向量的数据类型来进行选择。

我们以一个简单的求和函数作为例子

```r
double sumC(NumericVector v){

    double total = 0;
    int num = v.size();
    for (int i = 0; i< num; i++){
        total += v[i];
    }
    return total;
}
```

`sumC`这个函数读取一个数值向量，对其进行结合，返回一个浮点型标量。当然这个浮点型标量在R语言中就表现为一个长度为1的向量。

由于R语言数据类型和`C++`的不是一一对应，因此在处理过程中，有些时候需要对数据进行转换

```bash
std::string(x[i])
```

## 变量创建

我们可以通过下面这些方式创建向量

```c++
// v <- rep(0, 3)
NumericVector v (3);
// v <- rep(1, 3)
NumericVector v (3,1);
// v <- c(1,2,3) 
// C++11 Initializer list
NumericVector v = {1,2,3}; 
// v <- c(1,2,3)
NumericVector v = NumericVector::create(1,2,3);
// v <- c(x=1, y=2, z=3)
NumericVector v = NumericVector::create(Named("x",1), Named("y")=2 , _["z"]=3);
```

创建矩阵的方法如下

```C++
// m <- matrix(0, nrow=2, ncol=2)
NumericMatrix m1( 2 );
// m <- matrix(0, nrow=2, ncol=3)
NumericMatrix m2( 2 , 3 );
// m <- matrix(v, nrow=2, ncol=3)
NumericMatrix m3( 2 , 3 , v.begin() );
```

数据框和列表的创建依赖于已有向量，和R语言中创建形式相似。

```C++
// Creating DataFrame df from Vector v1, v2
DataFrame df = DataFrame::create(v1, v2);
// When giving names to columns
DataFrame df = DataFrame::create( Named("V1") = v1 , _["V2"] = v2 );

// Create list L from vector v1, v2
List L = List::create(v1, v2);
// When giving names to elements
List L = List::create(Named("name1") = v1 , _["name2"] = v2);a
```

## 元素获取和赋值

在选取元素之前一定要注意，`C++`是以0为基，而R是以1为基。

我们可以用`[]`和`()`进行数据选取向量中的元素并复制，支持利用数值向量或者逻辑向量来选取多个数据

```C++
// 新建向量
NumericVector v  {10,20,30,40,50};
// 设置向两名
v.names() = CharacterVector({"A","B","C","D","E"});
// 准备用于获取数据的向量
NumericVector   numeric = {1,3};
IntegerVector   integer = {1,3};
CharacterVector character = {"B","D"};
LogicalVector   logical = {false, true, false, true, false};
// 获取数据
double x1 = v[0];
double x2 = v["A"];
NumericVector res1 = v[numeric];
NumericVector res2 = v[integer];
NumericVector res3 = v[character];
NumericVector res4 = v[logical];
// 赋值
v[0]   = 100;
v["A"] = 100;
NumericVector v2 {100,200};
v[numeric]   = v2;
v[integer]   = v2;
v[character] = v2;
v[logical]   = v2;
```

对于矩阵而言，建议只用`()`进行数据选取，因为没有`[row_index, col_index]`的操作。

```C++
// 创建一个5x5的矩阵
NumericMatrix m( 5, 5 );
// 获取0,2的元素
double x = m( 0 , 2 );
// 选择第1行
NumericVector v = m( 0 , _ );
// 选择第3列
NumericVector v = m( _ , 2 );
// 选择第1-2行，3-4列
NumericMatrix m2 = m( Range(0,1) , Range(2,3) );
```

对于DataFrame和List，两则都只支持`[]`选择其中向量。

## 成员函数

成员函数(member function)是`C++`面向对象编程中的一个概念，通过调用成员函数可以对类进行操作。

举个例子，对于R语言而言，向量是可以有名字的，例如

```r
x <- c(a = 1, b = 2, c = 3)	
```

R语言是一个面向对象的编程语言，所以这种命名向量其实认为是一个向量拥有了一个名为`names`的属性而已。

在R语言中`attr`函数能够增加对象的属性，因此上面代码等价于

```r
y <- c(1,2,3)
attr(y, 'names') <- c('a','b','c')
```

而`C++`中的向量和矩阵拥有一个成员函数`.name()`就可以获取/修改/增加变量命名。

```C++
NumericVector addname2(NumericVector x, CharacterVector y){
    x.names() = y;
    return x;
}
```

更通用的是`.attr`函数,它能够修改和增加任意属性

不同数据类型拥有不同的成员函数，参考

- 向量:  https://teuder.github.io/rcpp4everyone_en/080_vector.html#member-functions
- 矩阵: https://teuder.github.io/rcpp4everyone_en/100_matrix.html#member-functions-1
- 数据框:  https://teuder.github.io/rcpp4everyone_en/140_dataframe.html#member-functions-2
- 列表: https://teuder.github.io/rcpp4everyone_en/150_list.html#member-functions-3
- 字符串:  https://teuder.github.io/rcpp4everyone_en/170_string.html#member-functions-4
- 日期: https://teuder.github.io/rcpp4everyone_en/180_date.html#member-functions-5
- 时间: https://teuder.github.io/rcpp4everyone_en/190_datetime.html#member-functions-6

## 缺失值

> 同R一样，缺失值之间以及缺失值和正常值之间的比较没有意义，因此在写代码过程时要提前考虑到缺失值的可能，或者是用R语言处理完缺失值，把不含缺失值的数据作为`C++`函数的输入。

R语言采用比特模式对每一种数据类型进行标注，也就是针对每一种数据类型（Integer，Numeric，Character，Logical）保留一个比特作为缺失数据的标签值。

为了在C++中处理R的缺失值，Rcpp提供了四个变量对应R的四种数据类型

- NA_INTEGER:  整型，Integer
- NA_STRING: 字符串，Character
- NA_LOGICAL：逻辑性，Logical
- NA_REAL：双精度浮点型， Numeric

写一个函数进行讲解下缺失值的使用。

```C++
double sumC(NumericVector v, bool na_rm){

    double total = 0;
    for (int i = 0; i< v.size(); i++){
        if( ! NumericVector::is_na(v[i])){
            total += double(v[i]);
        } else if (NumericVector::is_na(v[i]) & na_rm){
            continue ;
        } else {
            return NA_REAL;
        }
    }
    return total;
}
```

`NumericVector::is_na()`用于判断是否为缺失值，不同的数据类型有不同的`is_na`。如果`na_rm=FALSE`，那么返回缺失值，也就是`NA_REAL`.



## 参考资料

- [CharacterVector元素转string](https://stackoverflow.com/questions/23282849/rcpp-charactervector-why-arent-individual-elements-stdstrings)
- https://teuder.github.io/rcpp4everyone_en
- https://adv-r.hadley.nz/rcpp.html#rcpp-intro