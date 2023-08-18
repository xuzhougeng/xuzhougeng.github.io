---
title: 如何在服务器上安装最新的R
date: 2019-09-23 11:52:13.832
updated: 2020-04-13 12:51:50.085
url: /archives/install-r-in-linux
categories: Linux
tags: 软件安装
---



R语言在服务器上安装是一个比较可麻烦可简单的事情，这里记录下R语言在两个比较常见的Linux发行版的安装方法，分别是CentOS和Ubuntu。

## 通用方法（无需Root）

只要你的服务器能够安装conda，那么你就可以用conda去安装你的R语言。conda已经不再局限于最早的Python的环境管理了，而是扩展到R, Java, C/C++等编程语言。

> Package, dependency and environment management for any language—Python, R, Ruby, Lua, Scala, Java, JavaScript, C/ C++, FORTRAN, and more.

让我们使用`conda search r-base`在conda的频道中检索，

```bash
...
r-base                         3.5.1   hfb2a302_1009  anaconda/cloud/conda-forge
r-base                         3.5.1   hfb2a302_1010  anaconda/cloud/conda-forge
r-base                         3.6.0      hce969dd_0  pkgs/r
r-base                         3.6.1      h6e652e1_3  anaconda/cloud/conda-forge
r-base                         3.6.1      h8900bf8_0  anaconda/cloud/conda-forge
r-base                         3.6.1      h8900bf8_1  anaconda/cloud/conda-forge
r-base                         3.6.1      h8900bf8_2  anaconda/cloud/conda-forge
r-base                         3.6.1      hba50c9b_4  anaconda/cloud/conda-forge
r-base                         3.6.1      hce969dd_0  pkgs/r
```

发现能找到好几个版本的R语言，我推荐通过新建环境的方式安装不同版本的R语言，这样就能在不同环境间切换。

```bash
conda create -n r-3.6.1 r-base=3.6.1
```

之后用`conda activate r-3.6.1`调用R的环境即可。

这个方法的优点是不需要root权限，安装方便，不过听过在使用的时候或许会出现一些bug，我还没有遇到。

除了conda外，我们可以通过手工解决R语言的依赖环境，通过源码安装最新的R语言，这个方法也不依赖于平台。不过你需要看下这篇[无root权限解决编译时的依赖问题](/archives/Compile-Software-in-Linux-Without-Root), 但是很麻烦。还是建议用conda比较合适。

假如需要用Jupyter notebook调用R，那么安装方式为

```bash
conda install -c r r-irkernel 
# or
conda install -c r r-essentials
```

因为 conda install -c r r=3.6.x/r-base 默认不会安装 irkernel，而且先安装的 r=3.6.x/r-base 可能与后安装的 r-irkernel/r-essentials 产生冲突。

## CentOS

CentOS/RedHat是可以通过`sudo yum install R`的方式安装R语言，解决一切依赖问题，并且安装比较新的R版本

```bash
================================================================
 Package      Arch       Version     Repository    Size
================================================================
Installing:
 R            x86_64     3.6.0-1.el7  epel          30 k
Installing for dependencies:
 R-core       x86_64     3.6.0-1.el7  epel          57 M
 R-core-devel x86_64     3.6.0-1.el7  epel          109 k
 R-devel      x86_64     3.6.0-1.el7  epel           30 k
 R-java       x86_64     3.6.0-1.el7  epel           31 k
 R-java-devel x86_64     3.6.0-1.el7  epel           30 K
 ....
 Transaction Summary
=================================================================
Install  1 Package  (+373 Dependent packages)
Upgrade             (  14 Dependent packages)
```

这个方法安装的是比较新的R，基本上所有新版本能装的R包，它也能装了。但是如果你需要安装最新的R，那么就需要从头编译

先运行`sudo yum install R`搞定大部分依赖问题，然后你可能还得手动解决几个依赖问题，我遇到的是X11和libcurl

```bash
# X11
yum install xorg-x11-server-devel libX11-devel libXt-devel
# libcurl
yum install libcurl-devel
```

之后下载源代码编译安装

```bash
wget https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/base/R-3/R-3.6.1.tar.gz
tar -zxvf R-3.6.1.tar.gz
cd R-3.6.1
# --enable-R-shlib for Rstudio server
mkdir -p /opt/sysoft
./configure --enable-R-shlib --prefix=/opt/sysoft/R-3.6.1
make -j 8
make install 
```

其中`--enable-R-shlib`是Rsutdio-Server安装需要，而`--prefix`是指定安装路径

## Ubuntu

在Ubuntu上直接通过`sudo apt-get install r-base`的方式安装的不是最新版本的R，而是R-3.4版本。

> R 3.4 packages for Ubuntu on i386 and amd64 are available for all stable Desktop releases of Ubuntu prior to Bionic Beaver (18.04) until their official end of life date.

不过我们可以通过增加安装源的方式，使得能够通过`apt-get`的方式安装最新的R。

第一步，确认你的Ubuntu版本，是Xenial Xerus(16.04; LTS), Trusty Tahr (14.04; LTS), Bionic Beaver (18.04;LTS), Cosmic Cuttlefish (18.10), Disco Dingo (19.04)的哪一种。

第二步，根据你的服务器Ubuntu版本，按照需求复制下面的**其中**一行代码（一定要注意，是一行，不是全部复制）

```bash
deb https://cloud.r-project.org/bin/linux/ubuntu disco-cran35/
deb https://cloud.r-project.org/bin/linux/ubuntu cosmic-cran35/
deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/
deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/
deb https://cloud.r-project.org/bin/linux/ubuntu trusty-cran35/
```

然后用vim编译`/etc/apt/sources.list`, 添加你复制的内容到最后一行，我的服务器是 xenial，所以增加的是

```bash
deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/
```

除了修改`/etc/apt/sources.list`，还需要增加APT

```bash
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
```

之后，用下面的命令就可以安装最新的R

```bash
sudo apt-get update
sudo apt-get install r-base
sudo apt-get install r-base-dev
```

这个方法稍微麻烦些，据说通过这样子安装的R存在一些bug，不过我没有遇到。


假如你需要安装不同版本的R语言，那么就需要下载源代码进行编译安装。根据我的经验，你至少先得用下面这些命令安装R的依赖环境（可能还不够）

```bash
# 设置环境变量
export CFLAGS=" -fPIC" CXXFLAGS=" -fPIC"
## build-essential
sudo apt-get install -y build-essential
## java
sudo apt install -y openjdk-9-jdk
## 各种包
sudo apt install -y autoconf libz-dev libbz2-dev liblzma-dev libssl-dev 
# solve  libcurl problem
#sudo apt install -y libcurl4-openssl-dev  # not works for Ubuntu 16.04
sudo apt install -y libcurl4-gnutls-dev 
### curses
sudo apt-get install -y libncurses5-dev
### solve X11 problem
sudo apt-get install -y xorg-dev
### zlib2
wget http://zlib.net/zlib-1.2.11.tar.gz
tar -zxvf zlib-1.2.11.tar.gz && cd zlib-1.2.11 && ./configure && make && sudo make install && cd .. && rm -rf zlib-1.2.11
### bzip2
wget https://fossies.org/linux/misc/bzip2-1.0.8.tar.gz
tar -zxvf bzip2-1.0.8.tar.gz && cd bzip2-1.0.8  
# add -fPIC
sed -i 's/CFLAGS=/CFLAGS=-fPIC /' Makefile
make && sudo make install && cd .. && rm -rf  bzip2-1.0.8
```

> 假如你使用的是conda用户，那么安装之前，你需要用先退出conda环境，不然后续libcurl可能会一直提示出错

下载R的源代码，进行编译安装

```bash
wget https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/base/R-3/R-3.6.1.tar.gz
tar -zxvf R-3.6.1.tar.gz
cd R-3.6.1
# --enable-R-shlib for Rstudio server
./configure --enable-R-shlib --prefix=/opt/sysoft/R-3.6.1
make -j 8
make install 
```

## 其他

此外，我还推荐大家安装OpenSSL/Curl/XML/HD5F，可以节约之后的一点时间

Fedora, CentOS, RHEL

```bash
yum install openssl-devel
yum install libcurl-devel
yum install libxml2-devel
yum install hdf5-devel
```

Debian, Ubuntu等

```bash
sudo apt install libssl-dev 
sudo apt install libcurl4-openssl-dev
sudo apt install libxml2-dev
sudo apt-get install libhdf5-dev
```

