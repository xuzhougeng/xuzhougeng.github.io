---
title: 服务器上R调用png显示x11报错怎么办？
date: 2022-06-15 04:30:13.88
updated: 2022-06-15 04:55:56.555
url: /archives/how-to-solve-x11-error-in-server
categories: R
tags: 编程语言
---


太长不读版

- 治本的方法，服务器安装pango, 之后重新编译R语言
- 治标的方法，在R的配置文件中增加`options(bitmapType='cairo')`

服务器上装完R语言之后，发现自己的PNG函数无法正常调用，始终会出现一个X11连接失败的错误。

```r
> png('ab.png')
Error in .External2(C_X11, paste0("png::", filename), g$width, g$height,  :
  unable to start device PNG
In addition: Warning message:
In png("ab.png") : unable to open connection to X11 display ''
```

使用 `capabilities`查询之后，的确是没有启用X11。

```r
> capabilities()
       jpeg         png        tiff       tcltk         X11        aqua
       TRUE        TRUE        TRUE        TRUE       FALSE       FALSE
   http/ftp     sockets      libxml        fifo      cledit       iconv
       TRUE        TRUE       FALSE        TRUE        TRUE        TRUE
        NLS       Rprof     profmem       cairo         ICU long.double
       TRUE        TRUE       FALSE        TRUE        TRUE        TRUE
    libcurl
       TRUE
```

回溯到我源代码的配置过程，发现配置项目包括X11的界面支持，`Interfaces supported:X11, tcltk`. 但是我思索了下，这个X11支持指的应该是Linux的图形界面下的R，而不是终端的R，我的直觉告诉我，这个信息应该是和我们的报错无关了。

之后，我又去查了下其他资料，谢益辉提到这可能是png函数的type参数选择出现了问题，也就是下面这段函数

```r
    if (missing(type))
        type <- getOption("bitmapType")
```

进一步这个bitmapType选项的设置又和 `grDevices:::.onLoad`函数中的`.Call(C_cairoProps, 2L)`有关。其中.Call是R语言调用C语言编写函数的一种方式，回溯到对应的C代码(代码位置为src/library/grDevices/src/init.c)

```c
#ifndef _WIN32
/* This really belongs with the X11 module, but it is about devices */
static SEXP cairoProps(SEXP in)
{
    int which = asInteger(in);
    if(which == 1)
        return ScalarLogical(
#ifdef HAVE_WORKING_CAIRO
            1
#else
            0
#endif
            );
    else if(which == 2)
        return ScalarLogical(
#ifdef HAVE_PANGOCAIRO
            1
#else
            0
#endif
            );
    return R_NilValue;
}
#endif
```

还好，我略懂一些C语言，发现这个输出结果又和变量 HAVE_PANGOCAIRO 有关。 搜索R源代码中所有和HAVE_PANGOCAIRO有关的内容，发现在configure有相关记录，对应的逻辑判断语句如下

```text
if test "x${r_cv_has_pangocairo}" = xyes; then

printf "%s\n" "#define HAVE_PANGOCAIRO 1" >>confdefs.h

fi
```

进一步，我们在configure中查找 r_cv_has_pangocairo, 找到如下的判断语句

```text
else $as_nop
  if "${PKG_CONFIG}" --exists pangocairo; then
         r_cv_has_pangocairo="yes"
       else
         r_cv_has_pangocairo="no"
       fi

fi
```

那看来有可能是pangocairo这个库出现了问题。使用下列代码在我的服务器上测试，发现输出是no, 说明的确没有安装。

```bash
/bin/pkg-config --exists pangocairo && echo "yes" || echo "no"
```

我尝试着用下述代码安装pango

```bash
# CentOS
sudo yum install pango-devel pango libpng-devel
# ubuntu
sudo apt install libpango-1.0-0
```

对R进行重新编译安装，发现 `getOption("bitmapType")` 得到的结果是cairo，符合预期了。

当然，你可能也不需要那么麻烦的方法，还有一种简单的思路是在 `~/.Rprofile`中加入这一行，手动指定cairo。

```
options(bitmapType='cairo')
```

注意: 直接用apt/yum安装的R， `getOption("bitmapType")` 输出直接是 cairo.

参考资料

- https://stackoverflow.com/questions/17243648/cant-display-png
