---
title: 如何在不使用root权限下安装Latex
date: 2020-10-22 02:20:57.903
updated: 2020-10-22 02:20:57.903
url: /archives/install-latex-without-root
categories: 其他
tags: 软件安装
---

LaTex的安装并不需要管理员权限

> You do not need to be root (administrator on Windows) to install, use, or manage TeX Live. Indeed, we recommend installing it as a normal user, except perhaps on MacOSX, where it's conventional to install as administrator.

下面是Ubuntu的安装方法

第一步: 下载ISO文件(清华镜像源)

下载地址: [https://mirrors.tuna.tsinghua.edu.cn/CTAN/systems/texlive/Images/texlive2020-20200406.iso](https://mirrors.tuna.tsinghua.edu.cn/CTAN/systems/texlive/Images/texlive2020-20200406.iso)

> 其他镜像站点，可以在 [https://ctan.org/mirrors](https://ctan.org/mirrors) 进行选择

第二步: 在Windows上解压缩ISO文件，然后上传到服务器, 我们假设上传到家目录下的texlive 

第三步: 运行安装命令

```Bash
cd texlive
perl install-t 
```

输出信息如下

```Bash
======================> TeX Live installation procedure <=====================

======>   Letters/digits in <angle brackets> indicate   <=======
======>   menu items for actions or customizations      <=======

 Detected platform: GNU/Linux on x86_64
 
 <B> set binary platforms: 1 out of 6

 <S> set installation scheme: scheme-full

 <C> set installation collections:
     40 collections out of 41, disk space required: 6516 MB

 <D> set directories:
   TEXDIR (the main TeX directory):
     !! default location: /usr/local/texlive/2020
     !! is not writable or not allowed, please select a different one!
   TEXMFLOCAL (directory for site-wide local files):
     /usr/local/texlive/texmf-local
   TEXMFSYSVAR (directory for variable and automatically generated data):
     /usr/local/texlive/2020/texmf-var
   TEXMFSYSCONFIG (directory for local config):
     /usr/local/texlive/2020/texmf-config
   TEXMFVAR (personal directory for variable and automatically generated data):
     ~/.texlive2020/texmf-var
   TEXMFCONFIG (personal directory for local config):
     ~/.texlive2020/texmf-config
   TEXMFHOME (directory for user-specific files):
     ~/texmf

 <O> options:
   [ ] use letter size instead of A4 by default
   [X] allow execution of restricted list of programs via \write18
   [X] create all format files
   [X] install macro/font doc tree
   [X] install macro/font source tree
   [ ] create symlinks to standard directories
   [X] after install, set CTAN as source for package updates

 <V> set up for portable installation

Actions:
 <I> start installation to hard disk
 <P> save installation profile to 'texlive.profile' and exit
 <H> help
 <Q> quit

Enter command: 

```

我们需要输入D进入更改界面，然后输入1，输入我们目标安装路径，输入R回到主目录，最后输入I进行安装。

最后安装结果之后会有如下提示

```Bash
Add /home/xzg/opt/texlive/texmf-dist/doc/man to MANPATH.
Add /home/xzg/opt/texlive/texmf-dist/doc/info to INFOPATH.
Most importantly, add /home/xzg/opt/texlive/bin/x86_64-linux
to your PATH for current and future sessions.
Logfile: /home/xzg/opt/texlive/install-tl.log

```

即在 `~/.bashrc` 或 `~/.zshrc` 这些配置文件下添加如下内容

```Bash
export MANPATH=$MANPATH:/home/xzg/opt/texlive/texmf-dist/doc/man
export INFOPATH=$INFOPATH:/home/xzg/opt/texlive/texmf-dist/doc/info
export PATH=$PATH:/home/xzg/opt/texlive/bin/x86_64-linux
```

