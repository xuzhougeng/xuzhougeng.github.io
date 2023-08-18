---
title: 逃离dplyr:不使用group_by和arrange实现分组排序
date: 2020-01-09 16:07:03.766
updated: 2020-01-09 16:08:50.631
url: /archives/sort-by-group-without-dplyr
categories: R
tags: 
---

今天在写代码的时候，发现项目中出现了一些重复的代码，所以要把他们封装成一个单个函数。在封装的过程中，我遇到了一个让我头疼的问题。

在使用dplyr的时候，你可能会注意到一个非常有趣的细节，那就是你不用`""`来区别变量和字符串，`dplyr`能够帮你好这个事情。举个例子，下面的代码都是让`iris`数据集按照"Sepal.Length"进行排序。

```r
group_by(iris, Sepal.Length)
group_by(iris, "Sepal.Length")
```

这时候，让我们思考一个问题，如果在之前命名了一个`group.by <- "Sepal.Length"`，那么运行`group_by(iris, group.by)`的时候， 这个group.by会被替换成Sepal.Length吗？下面的代码会报错吗？大家可以思考一下，然后往下看。

```r
group.by <- "Sepal.Length"
group_by(iris, group.by)
```

实际上，运行上面的代码，你会得到一个报错"Error: Column `group.by` is unknown". `group_by`没有替换掉你的变量名。

为什么会出现这个情况？这个就涉及到dplyr编程的内容，具体可以参考[https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html](https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html)

参考dplyr的教程，如果要让上面的代码能够运行，我们需要需要在变量名前加上`!!`或者调用`UQ`函数

```r
group_by(iris, !!group.by)
group_by(iris, UQ(group.by))
```

假如你要创造一个函数，要用到`group_by`，那么你应该怎么写呢？我们的直觉就是下面的代码

```r
my_group_by <- function(df, group.by){
  df <- group_by(df,group.by )
}
```

根据前面的铺垫，你应该知道，运行`my_group_by(iris, Sepal.Length)`会出现报错，报错信息为" Error: Column `group.by` is unknown"。 于是你试着之前的解决方法加上了"UQ" 

```r
my_group_by <- function(df, group.by){
  df <- group_by(df, UQ(group.by) )
  return(df)
}
```

思考下，`my_group_by(iris, Sepal.Length)`能够得到结果吗？

很遗憾，代码报错了

```
 Error in splice(dot_call(capture_dots, frame_env = frame_env, named = named,  : 
  object 'Sepal.Length' not found 
```

正确的调用方法是`my_group_by(iris, "Sepal.Length")`.  当然由于你用习惯了dplyr，你希望是`my_group_by(iris, Sepal.Length)`调用代码，那么你的函数需要怎么写呢？为了解决这个问题，你可能要仔细阅读[https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html](https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html)才行，理解什么叫做"Quasiquotation"

经过一波折腾之后，我受够了这种将dplyr代码改成函数时的不一致性（也可能是我不熟练），我决定还是只用R基础代码来实现我的功能，不就是分组排序吗，为啥一定要用`group_by` 和 `arrange`呢.

无非就是先利用因子将数据库分成多个列表(split)，然后对每个列表按照某一列进行排序(lapply)，而这里排序过程就是获取从最大到最小的索引(order)，最后按行进行合并(do.call, rbind)而已呀。如下是实现的代码

```r
my.group.by <- function(df, group.by, sort.by,
                        decreasing = TRUE){
  
  df.split <- split(df, df[[group.by]])
  
  df.split.sort <- lapply(df.split, function(x){
    x.order <- order(x[[sort.by]],decreasing = decreasing)
    x <- x[x.order,]
    x
  })
   df <- do.call(rbind, df.split.sort)
  return(df)
}

my.group.by(iris, group.by = "Species", sort.by = "Sepal.Length")
````




