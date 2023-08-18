---
title: 服务器上安装RStudio-server
date: 2019-10-14 16:00:45.418
updated: 2019-10-14 16:00:45.418
url: /archives/Install-RStudio-server-in-Server
categories: Linux
tags: 服务器
---

# 服务器上安装RStudio-server

如果想在服务器上安装一个RStudio-server，你需要先保证自己拥有管理员权限，之后参考[如何在服务器上安装最新的R](/archives/Install-R-in-Linux)安装R语言，一定要**注意**在`./configure`的时候加上`--enable-R-shlib`，否则后续会出错。

RStudio-server分为两种版本，一种是开源免费版，另一个是商业专业版本。个人觉得两者最大的区别在于，商业版支持在多个版本的R语言之间进行切换，而开源免费版不行。

## CentOS篇

如果服务器安装的是CentOS/RedHat，那么需要保证它们的发行版本不等于6

从官方上下载rpm文件

```bash
wget https://download2.rstudio.org/server/centos6/x86_64/rstudio-server-rhel-1.2.5001-x86_64.rpm
```

如果是第一次安装，那么就是运行如下命令

```bash
sudo yum install rstudio-server-rhel-1.2.5001-x86_64.rpm
```

假如是需要升级RStudio-server，比如我原先的是`1.1.456`最新的是`1.2.5001`, 需要先暂停当前的服务并卸载

```bash
/usr/sbin/rstudio-server stop
yum remove rstudio-server
```

之后才是安装

```bash
sudo yum install rstudio-server-rhel-1.2.5001-x86_64.rpm
```

安装完成之后，我们可以通过修改`/etc/rstudio/rserver.conf`更改端口和R所在路径

```bash
www-port=8080 # 端口, 默认8787
www-address=0.0.0.0
rsession-which-r=/opt/sysoft/R-3.6.1/bin/R # 安装R的路径
```

修改完成之后，用`rstudio-server restart`重启服务，没有任何信息就表示安装成功了。

当然你要是不放心，你还可以用`rstudio-server verify-installation`来验证下，如果没有任何输出信息就表示安装成功，假如出现下面这条信息，意味着你需要先用`rstudio-server stop`先暂停服务。

```bash
Server is running and must be stopped before running verify-installation
```

其实最直接的方法就是直接访问"IP:端口"，能够出现RStudio的登陆界面就意味着安装成功了。

## Ubuntu篇

我没有一台Ubuntu系统的服务器，只有一台Windows 10电脑有一个Linux子系统安装的是 Ubuntu 16.04.6 LTS。

下载Deb文件

```bash
sudo apt-get install gdebi-core
wget https://download2.rstudio.org/server/trusty/amd64/rstudio-server-1.2.5001-amd64.deb
```

如果不是第一次安装，需要是升级已有的RStudio-server，那么也需要先停用并卸载已有的RStudio-server

```bash
sudo /usr/sbin/rstudio-server stop
sudo apt-get remove rstudio-server 
```

然后安装

```bash
sudo gdebi rstudio-server-1.2.5001-amd64.deb
```

如果是在Windows的子系统下安装，会出现如下的警告，允许访问即可。

![Windows中警告](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/10/1571039719464-278ec32cdcb94ca1ad185f592228e4e1.png)

如果能够打开[http://127.0.0.1:8787](http://127.0.0.1:8787/), 就说明安装成功了。

如果想修改RStudio-server的端口和调用R版本，参考CentOS篇

## 参考资料

- https://rstudio.com/products/rstudio/download-server/redhat-centos/
- https://rstudio.com/products/rstudio/download-server/debian-ubuntu/
- https://support.rstudio.com/hc/en-us/articles/216079967-Upgrading-RStudio-Server