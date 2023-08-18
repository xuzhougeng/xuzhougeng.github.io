---
title: 在Windows上配置类Unix环境(MSYS2)
date: 2021-06-30 16:59:53.772
updated: 2021-06-30 17:08:10.174
url: /archives/configure-msys2-in-windows
categories: 其他
tags: 小技巧
---

MSYS2通过精心整合一些已有的工具(如gcc, qt)，使得在Windows上也能使用类unix环境，调用类unix开发的软件。

它的官方地址为: [https://www.msys2.org/](https://www.msys2.org/)

> Cygwin尝试提供一个POSIX兼容环境，使得原本能在类unix环境里运行的软件，不需要经过太多修改就能直接在Windows上使用。它的特点是大而全，而MSYS2基于Cygwin提供优秀环境，使用了优秀的包管理工具pacman，更倾向于尝试提供一个构建原生Windows软件的环境。

## 安装步骤

第一步:在官网上下载MSYS2的[安装包](https://github.com/msys2/msys2-installer/releases/download/2021-06-04/msys2-x86_64-20210604.exe)

e地址: https://github.com/msys2/msys2-installer/releases/download/2021-06-04/msys2-x86_64-20210604.exe

第二步: 运行MSYS2的安装程序，注意记录的你安装路径，我是修改了默认的安装路径为 `D:\MSYS2`

![Fig1](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/06/image-5e76f5d89c3c40f482e795683817c8c9.png)

安装成功之后，你的菜单栏里就有了一个MSYS2的启动方式

![Fig2](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/06/image-589dc32f69ab4894a011834272bdcdfc.png)

第三步: 更新pacman的包数据库

```bash
pacman -Syu
#-S 表示同步 Sync
#u: -u --sysupgrade
#y: -y --refresh  
```

![Fig3](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/06/image-f87c788fabb74e61869ef3133fed1cde.png)

安装结束后，会关闭窗口，之后重新打开, 继续更新剩余的基础包

```bash
pacman -Su
```

第四步: 安装mingw-w64 GCC编译工具套装

```bash
pacman -S --needed base-devel mingw-w64-x86_64-toolchain
#  gcc vim cmake
pacman -S gcc vim cmake
# git 版本控制
pacman -S git 
```

中间会有一些选项，全部enter默认即可。

## 其他

### 在CMD中调用MSYS2编译程序

为了能够在cmd中使用MSYS2编译的软件，我们需要在**环境变量**中增加相应的路径，否则会出现如下报错

![Fig4](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/06/image-f4d657b9c21549f4a301e14f5cf34b99.png)

我安装的软件在D盘，因此我添加的路径是 `D:\msys64\usr\bin`

![Fig4](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/06/image-9618ca6bb31546b88221a6e9979a0775.png)

![Fig5](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/06/image-246638b5263e4f27b40d25a20df12db1.png)

### Pacman包管理

Pacman的常见操作

安装

```bash
pacman -S <package_names|package_groups>
```

移除

```bash
pacman -R  <package_names|package_groups>
```

搜索

```bash
pacman -Ss <name_pattern>
```

- package_names: 包的名字
- package_groups: 包的所在组
- name_pattern: 包的可能命名
