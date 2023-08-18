---
title: MacOS原生Arm64架构的R解决依赖编译问题
date: 2022-02-21 09:51:46.786
updated: 2022-02-21 10:01:03.216
url: /archives/how-to-compile-r-package-in-arm64-macos
categories: R
tags: MacOS
---


目前(截止到2022年2月份), Arm64架构的R能够直接安装CRAN绝大部分的R包，但是Bioconductor的R包都需要编译（这也就为什么我推荐用Intel架构的R，可以避免遇到要编译R包里的坑）。

编译R包最怕遇到的问题，就是遇到R以外的底层库的依赖，比如说Rhtslib, 就要求编译htslib, 而htslib依赖lzma是MacOS默认没有安装，因此我们会因为缺少xz库的支持，导致htslib编译失败，从而导致Rhtslib安装失败，间接导致所有依赖于Rhtslib的R包无法使用。还有一些R包编译需要gsl库，就会因为找不到 `gsl-config` 而报错。

经过我不断踩坑后，我整理了解决方案。后续，如果你也打算在Arm64的Mac上安装原生的R，就按照我的思路来配置


首先，你需要安装homebrew, 国内用户推荐按照清华镜像里的方案来配置，https://mirrors.tuna.tsinghua.edu.cn/help/homebrew/ , 能够避免网络问题导致的坑。在海外的朋友，就按照官方站点, brew.sh 的命令来配置。

homebrew会安装在 `/opt/homebrew`下，该目录里面包括后续安装的软件(bin), 头文件(include)和动态链接库(lib). 

之后，我们需要为R配置相关环节，使得R能够调用到homebrew安装的环境。

我们需要打开R，然后运行如下命令

```r
file.edit(file.path(Sys.getenv("R_HOME"), "etc", "Makeconf"))
```

如此就能打开R的Makeconf文件，R就通过该文件里面的参数来进行R包和底层依赖的编译。我们搜索搜索 CFLAG, 添加  `-I/opt/homebrew/inlucde`， 搜索LDFLAGS, 添加 `-L/opt/homebrew/lib`, 效果如下

```conf
CFLAGS = -I/opt/homebrew/include -falign-functions=64 -Wall -g -O2 $(LTO)
LDFLAGS = -L/opt/R/arm64/lib -L/opt/homebrew/lib
```

重启R/RStudio之后，你就能编译安装Rhtslib（当然你需要先运行 `brew install xz`，把lzma.h 装好`)

```R
BiocManager::install("Rhtslib")
```

此外你可能还会遇到`gsl-config not found`报错，尽管你明明用`brew install gsl`安装了软件，也能够在命令行中成功调用 `gsl-config`，但是R就是视而不见。 这是因为R读取的变量，默认**就来自于**MacOS的环境变量来自于 `/etc/path`文件和 `/ect/path.d`目录里。 因此，你即便在.zshrc, 或者.bashrc里面加入了 homebrew 的路径，R也是不认的。

解决方案，打开R，然后运行如下命令,

```R
file.edit("~/.Rprofile")
```

添加如下内容，修改R启动时的环境变量

```
Sys.setenv("PATH"=paste0( "/opt/homebrew/bin:/opt/homebrew/sbin:",  Sys.getenv("PATH")))

```

重启R/Rstudio之后，R就能够用到homebrew安装的程序了。

或者，我们也可以下载到本地，用 `R CMD INSTALL`进行编译安装，以DirichletMultinomial为例

```
wget https://mirrors.tuna.tsinghua.edu.cn/bioconductor/packages/3.14/bioc/src/contrib/DirichletMultinomial_1.36.0.tar.gz
R CMD INSTALL DirichletMultinomial_1.36.0.tar.gz
```

总结一下:

1. 通过homebrew来安装依赖环境
2. 编译R的Makeconf使得R能够调用homebrew里的include和lib
3. 编译R的.Rprofile修改PATH，使得R能够用到hombrew的bin的程序

