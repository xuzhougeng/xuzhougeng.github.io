---
title: R语言文件读取的基础知识
date: 2019-12-09 15:47:47.018
updated: 2019-12-09 15:47:47.018
url: /archives/R-basic-Input-Output-function
categories: R
tags: 
---

在使用Python读取文件的时候，通常我需要用到`open()`打开一个文件，接着是`readline()`或`readlines()`读取文件内信息，最后用`close()`关闭文件。但是在用R语言读取数据时，我就直接用到比较高层的函数，例如`read.table`或者`data.table::fread`。R语言也有比较底层的输入输出(I/O)函数，了解这些知识，可以让我们更好的读取文件。

## 连接(Connection)

连接(connect)是R语言中用于多种I/O操作的基本机制。R提供了9种函数用于为不同输入建立连接，file/url/gzfile/bzfile/xzfile/unz/pipe/fifo/socketConnection，基本上从名字上就能看出这些函数的使用场景。

举个例子，如果想要爬取一个网页，我们可以用`url`构建目标站点的连接

```R
url_con <- url(description="http://xuzhougeng.top", open="r", encoding = "UTF-8")
class(url_con)
[1] "url"        "connection"
```

该`url_con`就可以作为`read.table/readLines`的文件参数。

上述的`description`就是指用于将要被建立的连接，可以是一个文件名或URL地址，例子中用的就是一个URL。也可以是其他类型输入。

如果需要在命令行中使用R脚本，让R读取管道的输入，就是"stdin"

```R
std_in <- file("stdin", "r")
```

如果要读取剪切板的内容，则是"clipboard"

```r
clip_text <- file("clipboard", "r")
```

其最终返回一个"connection"对象。

## "特殊"连接

这里的"特殊"连接指的是不怎么用的到建立连接的方式，可能在某些特殊的情况下使用。

### 场景一: 从控制台(console)里获取输入

我们可以用`stdin()`来从console中获取输入

```r
> a <- readLines(stdin(), n=2)
hello
world
> a
[1] "hello" "world"
```

这里的`n=2`表示读取两行的输入。该方法不能用于shell命令行读取管道的输入

与`stdin()`对应的还有`stdout()`，`stderr()`和`nulffile()`，这几个命令很少使用。

### 场景二: 从复制的文本中获取输入

使用`clipboard`有些时候不清楚自己到底会从剪切板里读取什么样的输入，而且Linux命令行未必有X11支持，所以更好的方法是能够将文本粘贴到脚本中，然后进行读取。

为了解决这一问题，我们需要借助于`textConnection`将字符串中转成连接。

```r
tmp <- "a   b   c   d
1   1   1   B
2   1   2   C
3   1   3   C
4   1   4   C
5   1   5   B
6   1   6   B
7   1   7   C
8   1   8   A
9   1   9   B
10  1   10  C
"
df <- read.table(textConnection(object = tmp),sep="\t", header = T)
```

和`textConnection`类似的，还有一个`rawConnection`用于构建`raw`对象的连接

## 读取/输出函数

构建的连接取决于模式可以进行读取和输出操作。R自带函数中和读取/输出相关的可以分成两类

文本模式: readLines, writeLines, cat, sink, scan, parse, read.dcf, dput, dump

二进制模式: readBin, readChar, writeBin, writeChar, load, save

## 连接管理

对于连接对象，除了读写以外，还有一些列的函数用于查看和管理连接。

- `open`: 打开连接
- `close()`: 关闭连接
- `seek()`: 在连接中进行跳转
- `isOpen`判断连接是否打开
- `isIncomplete`判断连接的读取最后一次是否完整， 写出时是否有未输出的内容
- `flush`刷新连接的输出流，保证内容都完全写出。用`close()`安全的关闭连接
- `showConnections`: 显示目前的所有连接
- `getConnection`: 获取指定连接
- `closeAllConnections`: 关闭所有连接

大部分函数都不好找到使用场景，除非你要编写某些数据类型的读取函数。