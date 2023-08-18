---
title: 如何提高Github的访问速度
date: 2021-04-06 09:09:19.16
updated: 2021-04-06 09:27:28.388
url: /archives/how-to-accelerate-access-to-github
categories: Linux
tags: 小技巧
---

最近总觉得Github的访问速度变慢了，导致我的工作效率也肉眼可见的降低了，主要体现在代码数量和质量的双重降低。

为了解决这一问题，我通过网络检索找到了一个非常好的工具，叫做dev-sidecar (https://github.com/docmirror/dev-sidecar), 这个工具的好处在于，图形化界面，完美解决了windows和MacOS上的问题，但是问题也出在图形化界面上，这意味着没有配置图形界面的Linux就用不了。

但是这难不倒我们，因为我看到了这个工具的参考部分。

![reference](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/04/image-faa3d53abb1a485eb99513971e70ed3e.png)

从中, 我们可以看到一个非常重要的地址，fastgit, 这是一个对于 GitHub.com 的镜像加速器。它的使用方法及其简单，也就是将原来的github.com替换成 hub.fastgit.org即可。

以 bioconda/bioconda-recipes为例，在使用原版的GitHub时，我们的下载速度基本上维持在2MiB/s，某些时候可能不到100Kb.

```bash
$ time git clone https://github.com/bioconda/bioconda-recipes
Cloning into 'bioconda-recipes'...
remote: Enumerating objects: 64, done.
remote: Counting objects: 100% (64/64), done.
remote: Compressing objects: 100% (50/50), done.
remote: Total 276716 (delta 30), reused 31 (delta 11), pack-reused 276652
Receiving objects: 100% (276716/276716), 303.49 MiB | 3.72 MiB/s, done.
Resolving deltas: 100% (152345/152345), done.
Checking out files: 100% (17906/17906), done.
git clone https://github.com/bioconda/bioconda-recipes  46.51s user 7.04s system 53% cpu 1:39.24 total
```

但是使用fastgit加速之后，我的下载速度直接飙升到10Mib/s以上，峰值可以达到30Mib/s. 

```bash
$ time git clone https://hub.fastgit.org//bioconda/bioconda-recipes
Cloning into 'bioconda-recipes'...
remote: Enumerating objects: 64, done.
remote: Counting objects: 100% (64/64), done.
remote: Compressing objects: 100% (50/50), done.
remote: Total 276716 (delta 30), reused 31 (delta 11), pack-reused 276652
Receiving objects: 100% (276716/276716), 303.23 MiB | 21.16 MiB/s, done.
Resolving deltas: 100% (152348/152348), done.
Checking out files: 100% (17906/17906), done.
git clone https://hub.fastgit.org//bioconda/bioconda-recipes  44.88s user 6.09s system 101% cpu 50.367 total
```

速度的提升可能是其次的，最重要的是原本因为网络问题的fatal error导致根本下载不了的repos，现在起码能保证能克隆到本地了，实现了从0到1的进步。

当然，如果觉得每次都需要替换URL太过麻烦，fastgit还支持直接修改Git的配置，即

```bash
git config --global url."https://hub.fastgit.org/".insteadOf "https://github.com/"
git config --global protocol.http.allow always
```

之后，原本需要从Github上克隆的资源都会被定向到fastgit上，不再需要手动进行修改了。

不过这一缺陷在于，部分时候镜像站点可能会出现不可用的情况，此时你从Github克隆时依旧会被改向到镜像站点，你会误以为原站点出现了问题。不过只要我们人人都去支持下这个项目，随着可用节点的增多，这一缺陷也不是问题了。
