---
title: 编译安装容器singularity
date: 2023-08-27 12:36:26
update:
categories:
tags:
---

Singularity的编译不是完全需要root权限，但是有root权限会更加简单，如果你的服务器上没有singularity，用的是ubuntu，可以把这篇文章发给管理员，让他帮你装一个。

编译的Singularity的环境是Ubuntu 20.04，先安装必要的依赖环境

```bash
sudo apt-get update && sudo apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools
```

然后你需要安装大于1.20的GO，否则无法编译新的singularity，这里用的GO版本是1.21。

```bash
wget https://go.dev/dl/go1.21.0.linux-amd64.tar.gz 
sudo tar -C /usr/local -xzf go1.21.0.linux-amd64.tar.gz 
```

添加到自己的环境变量中

```bash
echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
source ~/.bashrc
```

接着，我们需要从github上克隆一份singularity到本地，并切换到对应目录中

```
mkdir -p $GOPATH/src/github.com/sylabs
cd $GOPATH/src/github.com/sylabs
git clone https://github.com/sylabs/singularity.git
cd singularity
```

克隆之后，还需要下载这个库中的一些GO依赖库

```bash
git submodule update --init
go get -u -v github.com/golang/dep/cmd/dep
```

最后就是编译安装

```
./mconfig
make -C builddir
sudo make -C builddir install
```

最后会得到两个文件，一个是执行文件`/usr/local/etc/singularity`，一个是配置文件`/usr/local/etc/singularity/singularity.config`。

需要注意的是，由于安装过程中，需要从GitHub上下载数据，因此要求服务器能够顺畅的访问GitHub，否则会因为网络原因失败。