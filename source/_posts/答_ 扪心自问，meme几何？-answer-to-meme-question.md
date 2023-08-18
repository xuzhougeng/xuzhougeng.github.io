---
title: 答扪心自问，meme几何？
date: 2019-12-11 11:46:14.345
updated: 2019-12-11 11:46:59.892
url: /archives/answer-to-meme-question
categories: R
tags: 可视化
---

在2018-04-06，Y叔推送了一篇文章[扪心自问，meme几何？](https://mp.weixin.qq.com/s/BPY4dPFrrNrDKccMvax07w)。从一个不到146行的`meme.R`出发，提了5个问题。让检查我们的ggplot2的理解程度。

从我第一次接触ggplot2开始，至今差不多过去了2年多时间。Y叔推荐的「ggplot2: 数据分析和图形语法」和「R绘图系统」也被我翻得书页泛黄，加上近日又在学习Hadley写的[extending-ggplot2](https://cran.r-project.org/web/packages/ggplot2/vignettes/extending-ggplot2.html)，近日终于有所悟，尝试解答Y叔提出的几个问题。

## 1. `meme`的输出不是图？为什么能画出图？

这一题，Y叔在[R绘图系统](https://book.douban.com/subject/26792674/)的书评中给出了答案。

> 感谢译者送我的签名版，这是最全面介绍R绘图系统的书，没有之一。Base因为有大量的绘图函数还大量在使用，但做为新人学习，必须是grid系统！因为图是对象，可以操作，只有在需要渲染成静态图片的时候才产生图片。

更加详细一点，核心语句在于`imageGrob <- rasterGrob(x)`将读取的图像转成了grid的Grob对象，之后在此基础上构建了`p`，它继承了`meme`和`recordeplot`。这里的**继承**是关键词，也就是第二题的答案所在。

## 2. 为什么`+`可以改变图的内容和状态？

解答这一题需要两个关键知识

1. R语言是一个面向对象编程的语言，里面有一类泛型函数，可以根据你的对象类型自动调用对应的函数。
1. `+`是函数。

Y叔为`meme`对象专门定义了一个泛型函数`+.meme`, 因此在调用`+`的时候，也就是调用了`+.meme`函数。

## 3. 为什么`ggsave`能识别`meme`对象

这一题是讨论`ggsave`的本质，如果你直接在命令行里敲`ggsave`，他会输出`ggsave`的源代码，倒数第二句就是答案所在。

```r
grid.draw(plot)
```

`ggsave`在绘图商调用的是`grid.draw`， 这是用来绘制一个`grob`对象。而无论是cowplot，还是meme, 它们都建立在grid系统下，也就能够用`grid.draw`画出来。

而如果你调用的是`print(meme)`，那么泛型函数会尝试调用`print.grob`

## 4. 为什么使用传统的出图方式来画`meme`，在循环中需要显示`print(object)`？而`ggsave`则不用？到底区别在那里？

 这个问题稍微比较复杂， 我们需要先来实际的代码进行演示。

下面的for循环中，图形设备中不会出现图片，并且test.pdf打开的时候会显示图形损坏。

```r
library(meme)
u <- system.file("angry8.jpg", package="meme")

for ( i in 1:10){
  meme(u, "code", "all the things!")
}

# 图片输出到新的pdf中
pdf("test.pdf")
for ( i in 1:10){
  meme(u, "code", "all the things!")
}
dev.off()
```

而下面的代码中，图形设备会打印图片，并且test2.pdf能出现图片。

```r
for ( i in 1:10){
  p <- meme(u, "code", "all the things!")
  print(p)
}
# 图片输出到新的pdf中
pdf("test2.pdf")
for ( i in 1:10){
  p <- meme(u, "code", "all the things!")
  print(p)
}
dev.off()
```

为什么要在for循环里要用到`print`才行呢？

我们在R的控制台(console)运行meme时，实际上R会给你调用对应print函数答应。而在for循环中，它不会调用`print`。因此你必须要显示的调用`print`才行。

> In a loop, automatic printing is turned off, as it is inside a function

参考: https://stackoverflow.com/questions/4716152/why-do-r-objects-not-print-in-a-function-or-a-for-loop

## 5. 为什么`meme`对象能够被`ggimage`和`cowplot`识别？

Y叔说的马甲其实就是指`meme`继承了`recordedplot`，不过现在版本的`cowplot`似乎搞不定

```r
cowplot::plot_grid(p, p, ncol=1, labels = c("A", "B"))
# 报错如下
Error in value[[3L]](cond) : 
  invalid "recordedplot": Incompatible graphics state
In addition: Warning messages:
1: In restoreRecordedPlot(x, reloadPkgs) :
  snapshot recorded in different R version (pre 3.3.0)
2: In doTryCatch(return(expr), name, parentenv, handler) :
  snapshot recorded with different graphics engine version (pre 11 - this is version 12)
```

而为什么`ggimage`能够识别呢？`ggimage`是创建了一个专门的`geom_image`图层，为此Y叔利用ggproto，基于grid系统创造了一个`GeomImage` 类。这个图层就是用来绘制图片。