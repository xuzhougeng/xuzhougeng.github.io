---
title: 如何向GitHub上传超过100MB的大文件
date: 2021-11-18 06:30:09.324
updated: 2021-11-18 06:31:48.354
url: /archives/upload-large-file-to-github
categories: 生信软件工具箱
tags: GitHub
---


Git是一个版本控制工具，Github是一个项目托管网站。一般来说，我们会对代码进行版本控制，而记录代码的文件通常也很小，所以，我一般也不会用这种上传大文件的需求。并且，往GitHub上传超过50MB的文件时，它会警告你，认为这会影响性能，超过100Mb，直接报错。综上，我觉得往GitHub上传大文件是我不需要掌握的技能。

直到最近，我发现我上传的一个项目里面居然有一个400Mb的数据时，我就不得不去处理之前我觉得没有必要的蠢问题。（总不可能把数据放在百度网盘让人下载吧，这不太体面）

往Github上传大文件，我们需要使用GitHub的LFS(Large File Storage)服务，其中免费版允许2GB(超过2GB就需要付流量费）。

首先，点击Github仓库的Setting中，

![Setting](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/11/image-dad3669d568e43afaca5aa7b5805b497.png)

在其中的Archives中勾选Git LFS对应选项

![Git LFS](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/11/image-aa8a4535672a46c3b214b6de580dd880.png)

第二步，在服务器上安装Git LFS

```bash
# centos
sudo yum install git-lfs
# ubuntu
sudo apt install git-lfs
```

第三步: 设置LFS

```
git lfs install
# 提示信息如下
Git LFS initialized.
``` 

第四步: 在你项目下设置你需要跟踪的大文件，例如我需要跟踪R语言保存的二进制文件，一般是以Rds, Rds, Rdata, Rda结尾的文件

```bash
git lfs track "*.rds"
git lfs track "*.Rds"
git lfs track "*.Rdata"
git lfs track "*.Rda"
# 让git跟踪记录lfs信息的文件
git add .gitattributes
```

第五步：正常使用Git管理你的项目，上传到GitHub

```bash
git add data/*.rds
git commit -m "add rds"
git push origin master 
``

上传的时候，提示信息如下

```nbash
Uploading LFS objects:   0% (0/1), 3.4 MB | 129 KB/s
```

整体上，并不复杂，相当于在基础款的Git上安装一个插件而已。

参考资料:

- https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github
