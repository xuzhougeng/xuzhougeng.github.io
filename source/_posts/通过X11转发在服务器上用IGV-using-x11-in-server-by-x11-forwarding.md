---
title: 通过X11转发在服务器上用IGV
date: 2021-01-18 10:07:16.696
updated: 2021-01-18 10:07:16.696
url: /archives/using-x11-in-server-by-x11-forwarding
categories: 生信软件工具箱
tags: 软件安装
---

通常情况下，我们无法是直接在服务器上打开IGV，因为它既要求服务器上**安装X服务**，本地也要安装**X服务**，才能将服务器上的图形信息转发到本地。

> 事实上是要配置好X11转发，就可以在服务器上运行其他的图形界面程序。

第一步: 在服务器上配置Xserver(需要管理员权限)

```Bash
sudo yum install  xorg-x11-server-Xorg xorg-x11-xauth xorg-x11-apps -y
```

第二步: 在windows上配置X服务

我们从 [https://sourceforge.net/projects/xming/](https://sourceforge.net/projects/xming/) 下载Xming

![Xming](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/01/image-bbb12dad77a0414d83c71915da50655c.png)

安装完成后，右下角有Xming Server的地址信息，一般会是0.0或者1.0

![Xnubg Server](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/01/image-a4f1114d2222410e9be8277e971a1923.png)

第三步: 设置Xshell或Putty

!配置Xshell](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/01/image-8dbd12df088e4b688b71f4b72f86ecd2.png)

![配置Putty](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/01/image-9ee5bdc814074cac90183d424e0e6f55.png)

第四步: 下载IGV并解压缩

```Bash
wget https://data.broadinstitute.org/igv/projects/downloads/2.8/IGV_Linux_2.8.13_WithJava.zip
unzip IGV_Linux_2.8.13_WithJava.zip 
```

运行该目录下的 `igv_hidpi.sh ` 或 `igv.sh`即可