---
title: Bioconductor安装的网络问题
date: 2022-12-31 07:35:27.426
updated: 2022-12-31 07:35:27.426
url: /archives/the-network-problem-when-using-bioconductor
categories: R
tags: biocondutor
---

由于Bioconductor的网站在国外，因此一部分用户可能在使用BiocManager安装R包时因网络问题而失败。

第一个问题是，`BiocManager::install()`在运行时会下载 http://bioconductor.org/config.yaml 然后检查当前系统是不是符合需求，R包有没有过期。

尽管config.yaml文件非常小，但是依旧有一小部分人会因为访问不了biconductor.org, 导致这步花费非常久的时间。

两个解决方案：

第一个是设置BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS来避免自检

```r
options(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS=FALSE)
```

第二个是，下载好这个config.yaml到本地, 设置BIOCONDUCTOR_CONFIG_FILE为该文件在本地的路径

```r
options(BIOCONDUCTOR_CONFIG_FILE="/路径/到/你下载的/config.yaml")
```

第二个问题是，用官方源下载速度非常慢，所以可以使用国内镜像，如

```r
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
```

亦或者用 `chooseBioCmirror()` 选择。

