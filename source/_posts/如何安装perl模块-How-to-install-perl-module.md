---
title: 如何安装perl模块
date: 2019-08-23 12:57:37.759
updated: 2019-08-23 12:58:12.629
url: /archives/How-to-install-perl-module
categories: Linux
tags: Perl | 软件安装
---

由于生物信息早期最多用的语言是perl，因此不可避免就要用别人的perl脚本或者基于perl的项目来处理数据。

使用perl脚本和使用其他编程语言的脚本没啥不同，毕竟你只要传入参数，它就能给你结果。因此对于我们这些不用perl写脚本，只需要调用的人而言，唯一要学会的事情就是**如何安装perl的模块”。

关于perl模块安装，最古老的方法就是使用`perl -MCPAN -e shell`或者是`cpan`（两者等价），这也是我最先接触的方法，这里介绍如何使用`local::lib`和`cpanm`实现非root权限安装perl模块。

## 使用系统自带的perl

安装任何软件最怕遇到的问题就是权限问题，因此我们需要先安装`local::lib`，使得我们能够将perl模块安装到任何地方，简单的说就是安装到我们的家目录下

第一步，下载源代码进行编译安装

```bash
wget https://cpan.metacpan.org/authors/id/H/HA/HAARG/local-lib-2.000024.tar.gz
tar xf local-lib-2.000024.tar.gz
cd local-lib-2.000024
perl Makefile.PL --bootstrap=~/opt
make test && make install
```

第二步：**设置环境变量**，使得perl在安装模块的时候会优先使用我们指定的路径

```bash
echo 'eval "$(perl -I$HOME/opt/lib/perl5 -Mlocal::lib=$HOME/opt)"' >> ~/.bashrc
```

先用`perl -I$HOME/opt/lib/perl5 -Mlocal::lib=$HOME/opt`表示运行前先添加`$HOME/opt/lib/perl5`到自己的搜索路径`@INC`中，然后传入参数`$HOME/opt`执行模块`local::lob`，这个模块的执行结果会输出如下内容

```bash
Attempting to create directory /home6/wangjw/opt
PATH="/home/zgxu/opt/bin${PATH:+:${PATH}}"; export PATH;
PERL5LIB="/home/zgxu/opt/lib/perl5${PERL5LIB:+:${PERL5LIB}}"; export PERL5LIB;
PERL_LOCAL_LIB_ROOT="/home/zgxu/opt${PERL_LOCAL_LIB_ROOT:+:${PERL_LOCAL_LIB_ROOT}}"; export PERL_LOCAL_LIB_ROOT;
PERL_MB_OPT="--install_base \"/home/zgxu/opt\""; export PERL_MB_OPT;
PERL_MM_OPT="INSTALL_BASE=/home/zgxu/opt"; export PERL_MM_OPT;
```

这些就作为`eval`的参数进行执行，也就是说你重启终端后后，`PERL5LIB` `PERL_LOCAL_LIB_ROOT`,`PERL_MB_OPT`,`PERL_MM_OPT`这几个变量就会重新设置，以此保证你后续安装perl模块时，会优先安装到自己的选择的目录

第三步：安装cpam. 由于之前已经配置了`local::lib`，因此perl编译的工具都会默认安装到`~/opt`目录下

```bash
wget https://cpan.metacpan.org/authors/id/M/MI/MIYAGAWA/App-cpanminus-1.7043.tar.gz
tar xf App-cpanminus-1.7043.tar.gz
cd App-cpanminus-1.7043
perl Makefile.PL
make test && make install
```

第四步：使用国内镜像提高下载速度，可以通过别名的方式实现

```bash
echo 'alias cpanm="cpanm --mirror http://mirrors.163.com/cpan --mirror-only"' >>~/.bashrc
```

之后便可以使用`cpanm Module::Name`安装任意的软件了。

## 自己编译一个perl

自己编译Perl的好处就在于之后的perl模块都会安装到自己的Perl目录下，而不会对系统造成影响。

```bash
cd ~/src
wget -4 http://www.cpan.org/src/5.0/perl-5.26.1.tar.gz
tar xf perl-5.26.1.tar.gz
cd perl-5.26.1
./Configure -des -Dprefix=$HOME/opt/sysoft/perl-5.26.1
make test
make install
```

然后用`perl -e '{print "$_\n" foreach @INC}'`会发现perl只会在自己的目录`~/opt/sysoft/perl-5.26.1`下查找模块。那么使用`cpanm Module::Name`安装的任何包都只会安装到`~/opt/sysoft/perl-5.26.1`下，你也不需要安装`local::lib`了

## conda的perl和系统的perl冲突

有一次我遇到这个问题

```bash
perl: symbol lookup error: /home/wangjw/perl5/lib/perl5/x86_64-linux-thread-multi/auto/Cwd/Cwd.so: undefined symbol
```

这个问题是因为用系统perl安装的软件被conda的perl优先查找到导致，用`perl -V`和`perl -e '{print "$_\n" foreach @INC}'`可以发现conda的perl查找路径低于我为系统perl安装的路径，解决方案如下

```bash
export PERL5LIB=""
```

