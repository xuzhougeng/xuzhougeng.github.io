---
title: 如何在MacOS/Linux的RStudio上切换R版本
date: 2019-11-22 13:04:06.818
updated: 2019-11-22 13:04:06.818
url: /archives/Using-different-R-version-in-RStudion-linux-macOS
categories: R
tags: 
---

在MacOS或Linux上使用RStudio有一点不足，就是不能方便在RStudio上切换R版本。

虽然不方便，但还是可以切换。

首先，我们需要通过编译到方式安装不同版本的R，参考[如何在服务器上安装最新的R](/archives/Install-R-in-Linux)

我习惯将软件安装在`/opt/sysoft`下

```bash
sudo mkdir -p /opt/sysoft
sudo chown `whoami` /opt/sysoft
```

下载源代码，并安装

```bash
wget https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/base/R-3/R-3.5.3.tar.gz
tar xf R-3.5.3.tar.gz
cd R-3.5.3
./configure --enable-R-shlib --prefix=/opt/sysoft/R-3.5.3 --with-x=no
make -j 8
make install
#
wget https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/base/R-3/R-3.6.1.tar.gz
tar xf R-3.6.1.tar.gz
cd R-3.6.1
./configure --enable-R-shlib --prefix=/opt/sysoft/R-3.6.1 --with-x=no
make -j 8
make install
```

在安装完成之后，之后切换R语言，只需要用不同版本的R覆盖`/usr/local/bin`下面的R即可。

```bash
ln -sf /opt/sysoft/R-3.5.3/bin/R /usr/local/bin
ln -sf /opt/sysoft/R-3.6.1/bin/R /usr/local/bin
```

之后打开RStudio就是更换版本后的R


