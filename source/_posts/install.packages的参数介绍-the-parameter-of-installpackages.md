---
title: install.packages的参数介绍
date: 2020-04-11 16:22:07.266
updated: 2020-04-11 16:22:07.266
url: /archives/the-parameter-of-installpackages
categories: R
tags: 
---

今天解决解决了一个R包安装的问题，并且硬着头皮把`install.packages`和`download.file`的说明从头到位看了一遍，应该再也没有一个R包安装能为难到我了。

## 案例

### 问题描述

能够用浏览器访问镜像站点，但是在安装R包时遇到如下问题，

```bash
# CRAN
Warning in install.packages :
  unable to access index for repository https://mirrors.ustc.edu.cn/CRAN/src/contrib:
  cannot open URL 'https://mirrors.ustc.edu.cn/CRAN/src/contrib/PACKAGES'
Warning in install.packages :
  unable to access index for repository https://mirrors.ustc.edu.cn/CRAN/src/contrib:
  cannot open URL 'https://mirrors.ustc.edu.cn/CRAN/src/contrib/PACKAGES'
Warning in install.packages :
  package ‘ggtree’ is not available (for R version 3.5.1)
Warning in install.packages :
  unable to access index for repository https://mirrors.ustc.edu.cn/CRAN/bin/windows/contrib/3.5:
  cannot open URL 'https://mirrors.ustc.edu.cn/CRAN/bin/windows/contrib/3.5/PACKAGES'
# Bioconductor
Error: Bioconductor version cannot be validated; no internet connection?
In addition: Warning messages:
1: In library(package, lib.loc = lib.loc, character.only = TRUE, logical.return = TRUE,  :
  there is no package called ‘GSEABase’
2: In file(con, "r") : InternetOpenUrl failed: '??'
3: In file(con, "r") : InternetOpenUrl failed: 'on'
```

### 解决思路

第一步，确认R能否真的能够下载数据。检索到R用`download.file`进行文件下载，

```bash
 download.file(url = "https://upload-images.jianshu.io/upload_images/2013053-6e5c996e3a0d4c93.png",
              destfile = "test.png")
```

发现无法直接下载内容，证明R在连接网络时出现了问题

```bash
trying URL 'https://upload-images.jianshu.io/upload_images/2013053-6e5c996e3a0d4c93.png'
Error in download.file(url = "https://upload-images.jianshu.io/upload_images/2013053-6e5c996e3a0d4c93.png",  : 
  cannot open URL 'https://upload-images.jianshu.io/upload_images/2013053-6e5c996e3a0d4c93.png'
In addition: Warning message:
In download.file(url = "https://upload-images.jianshu.io/upload_images/2013053-6e5c996e3a0d4c93.png",  :
  InternetOpenUrl failed: 
```

第二步，根据报错信息, "InternetOpenUrl failed"进行检索，找到一种解决思路，也就是指定R访问网络的方法为`libcurl`

```bash
download.file(url = "https://upload-images.jianshu.io/upload_images/2013053-6e5c996e3a0d4c93.png",
              destfile = "test.png", methods="libcurl")
```

能够解决问题。

## 深入学习install.packages()

为了让自己能够更好解决R包安装问题，需要深入学习`install.packages`的原理(`BiocManager::install`本质上也是调用`install.packages`)。

先仔细阅读`install.packages()`的参数:

`pkgs`: 默认是包名，比如说"Matrix", 会自动从CRAN上检索对应的包，然后进行下载。如果你希望指定安装本地包，或者一个具体的网络地址，参考代码如下:

```R
# from url resource
install.packages("https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/windows/contrib/3.5/Matrix_1.2-15.zip", repos=NULL)
# from local
install.packages("~/../Desktop/Matrix_1.2-15.zip", repos = NULL)
```

- `lib`: R包放在那里，和`.libPaths()`有关

- `repos`: 镜像地址，当设置为NULL时，就可以安装本地包，或一个URLs。

- `contriburl`: 这个参数不常用，一般是你自己搞了一个本地镜像点时使用，该参数会覆盖掉`reops`. 与`type="both"`不兼容

- `method`: R包下载的方法，默认参数是`default`. 对于`file://`会调用`internal`，对于`ftps://`会调用`libcurl`,对于 `http://`,`https://`,`ftp://`, windows默认使用`wininet`,对于Unix类系统，默认使用`libcurl`. **注意**：如果Windows上用`capabilities("libcurl")`返回时TRUE, 那么也可以用`libcurl`,  Unix类系统无法使用`wininet`.

- `available`: 可以是`avaiable.packages`返回的镜像点中可用R包，也可以设置为NULL(这时函数内部会自动调用`avaiable.packages`). 与`type="both"`不兼容

- `destdir`: 下载的R包存放位置，`NULL`表示放在临时文件夹中，在关闭R后会被删除。

- `dependencies`:默认是NA，表示 `c("Depends", "Imports", "LinkingTo")`, TRUE表示对于要安装的R包是`c("Depends", "Imports", "LinkingTo", "Suggests")`依赖，依赖的依赖是`c("Depends", "Imports", "LinkingTo", "Suggests")`. **注意**： 对于二进制包，都会忽略"LinkingTo"

- `type`:  下载的是二进制包("binary")还是源代码"source". 如果设置为"binary", 依旧会先去检查该软件包最新的版本是否只有源代码，可用`options(install.packages.check.source = "no")`关闭。当设置为"source"时，只有不含"C/C++/Fortran"代码的R包可以被编译，如果R包中有`C/C++/Fortran`代码，那么Windows就需要安装Rtools。**注意**： 在Windows编译R包时，有一小部分需要设置`INSTALL_opts = "--force-biarch"` 或` INSTALL_opts = "--merge-multiarch" `, 建议后者。

- `configure.args`: 该参数只在源代码编译时使用，会传入`R CMD INSTALL`中

- `configure.vars`: 该参数只在源代码编译时使用， 类似于`configure.args`, 效果是在运行`configure`前设置环境变量。

- `clean`: 在`R CMD INSTALL`中加入`--clean`参数，用于清除临时中间文件。

- `Ncpus`： 编译时用多少CPU，加快编译速度。

- `verbose`: 是否输出安装时的信息。

- `libs_only`: 是否只安装64位或者32位的动态链接库

- `INSTALL_opts`: 源代码编译时`R CMD INSTALL`的额外传入参数，例如`c("--html", "--no-multiarch", "--no-test-load").`

- `quiet`: 安静模式，降低输出的信息量

- `keep_outputs`: 是否在当前工作目录下保留源代码编译后的输出文件。

- `...`的额外的参数来自于`download.file`, 主要就是`cache=TRUE`表示服务端缓存。默认是"TRUE"，如果是http://`和`https://`更建议用`cacheOK=FALSE`, 避免一些报错。

### 关于Secure URLs

对于`https://`或者`ftps://`这类URL，R在下载数据时会尝试对网站的证书进行验证。通常会调用操作系统安装的CA root certificates完成。

对于Windows系统，`method="libcurl"`时可能会出现问题，也就是Windows系统不提供有效的CA certificate bundle, 也就是说默认情况下，Windows的certificates是没有被验证过的。也就是`Sys.getenv("CURL_CA_BUNDLE")`返回结果为空，建议`Sys.setenv(CURL_CA_BUNDLE=file.path(Sys.getenv("R_HOME"),"/etc/curl-ca-bundle.crt"))`打开验证。

可以从<https://raw.githubusercontent.com/bagder/ca-bundle/master/ca-bundle.crt>下载curl-ca-bundle.crt的备份。

### 关于代理

 `wininet`调用系统中的`Internet Option`处理代理(proxy). 或者用`Sys.setenv`设置环境变量`http_proxy`,`ftp_proxy`

![互联网选项Internet Options](https://upload-images.jianshu.io/upload_images/2013053-8ed660914b813a4d.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

## 解决报错的一些建议

建议1：在用`install.packages`安装R包之前，先用`update.packages()`升级一下已有R包。

建议2：对于文章开头的报错，请修改`~/.Rprofile`，增加如下内容

```bash
options("download.file.method"="libcurl")
# # getOption("download.file.method")
options("url.method"="libcurl")
# getOption("url.method")
# options(internet.info = 0) # 进行HTTP传输的诊断，默认是2，只会显示最后的出错信息
```

也就是Windows尽量设置默认的method为libcurl，因为wininet未必一直支持HTTPS。 --<https://github.com/r-lib/remotes/issues/45#issuecomment-262955721>

建议3：如果遇到编译源代码报错，请在`install.packages()`中设置参数`INSTALL_opts = "--merge-multiarch"`和`clean=TRUE`

建议4：如果有些R包在安装时不能正确处理依赖关系，请在`install.packages()`中设置参数`depencies=TRUE`



