---
title: 在windows下编译bwa和samtools
date: 2021-06-30 17:07:37.498
updated: 2021-06-30 17:08:20.585
url: /archives/compile-bwa-and-samtools-in-windows
categories: 其他
tags: 小技巧
---


通过[在Windows上配置类Unix环境(MSYS2)](/archives/configure-msys2-in-windows)里的操作，我们成功的在Windows里拥有了一个类Unix的环境，后续，我们就可以在这里面编译一些软件，这里以bwa和samtools为例

以下是一些必须的依赖项

- zlib-devel
- git
- make
- gcc

用pacman安装依赖项，如zlib-devel和git

```bash
pacman -S zlib-devel
# 查漏补缺
pacman -S --needed base-devel 
pacman -S mingw-w64-x86_64-toolchain
pacman -S git 
pacman -S gcc vim cmake
```

如果连接git速度有问题，可以尝试按照如下方法配置镜像

```bash
git config --global url."https://hub.fastgit.org/".insteadOf "https://github.com/"
git config protocol.https.allow always
```

## 编译bwa

bwa只依赖于zlib，因此会很顺利

```bash
git clone https://github.com/lh3/bwa
# compile
cd bwa
make 
```



## 编译htslib

samtools依赖于htslib，因此需要先编译htslib

```bash
paman -S libbz2-devel liblzma-devel

git clone https://github.com/samtools/htslib
cd htslib
git submodule update --init --recursive 
# 不需要root
./configure
make -j 
make install
```

## 编译samtools

已知我不需要用到samtools 的 tview功能，则可以设置 `--without-curses `

```bash
git clone https://github.com/samtools/samtools
cd samtools
autoheader
autoconf -Wno-syntax  
./configure --without-curses 
make -j
make install 
```

如果需要直接在cmd中调用msys2编译的程序，需要将 MSYS2 安装目录下的 `usr/bin`添加到Windows的环境变量中Path，参考[在Windows上配置类Unix环境(MSYS2)](/archives/configure-msys2-in-windows)里的操作

可能报错:

```bash
make: *** 没有规则可制作目标“../htslib/htscodecs/htscodecs/arith_dynamic.c”，由“../htslib/hts-object-files” 需求。 停止。

```

原因: htslib没有克隆完全。

