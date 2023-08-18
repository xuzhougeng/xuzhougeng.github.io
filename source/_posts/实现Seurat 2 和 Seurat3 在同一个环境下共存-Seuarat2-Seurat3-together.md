---
title: 实现Seurat 2 和 Seurat3 在同一个环境下共存
date: 2019-08-14 10:31:49.412
updated: 2019-08-27 08:11:28.553
url: /archives/Seuarat2-Seurat3-together
categories: R
tags: Seurat
---

软件升级虽然是一件值得高兴的是，但是代码变化太大却不是一件好消息。比如说Seurat，这个单细胞分析最常用的R包，它的2.x版本和3.x版本的变化就是翻天覆地。

为了能够重现别人的代码，你可能需要重装2.3.4版本的Seurat，官方提供了安装脚本

```r
source("https://gist.githubusercontent.com/satijalab/beb9bb50dedc75ee023bd5d9be5fe684/raw/e103577735a2fba9da2ccca14ce1ac33e46c1bc4/old_seurat.R")
```

但是在window下安装会出一个报错, 提示无法在Windows下安装MacOS的版本

这个时候需要手动安装Seurat，注意，这里要求有Rtools才能进行编译

```r
install.packages('Seurat', repos = 'https://satijalab.org/ran', type = "source")
```

可以设置`install.packages`里的参数lib.loc，让Seurat安装到其他的文件加，就不会替换原来的Seurat，同时用`library`加载的时候，也需要设置`lib.loc`.

上面是简单可行易操作的方法，唯一的问题是你不能同时加载两个Seurat版本，当然你也不会想这样子做，所以这是一个伪需求。

下面就是瞎折腾环节， 我要将Seurat的2.3.4版本单独搞出一个R包，Seurat2，这样子就可以同时加载这两个R包。

下载Seurat 2.3.4

```bash
wget https://satijalab.org/ran/src/contrib/Seurat_2.3.4.tar.gz
```

解压缩它

```bash
tar xf Seurat_2.3.4.tar.gz
cd Seurat
```

删除MD5文件，因为它会做文件检验

```bash
rm MD5
```

修改里面所有的Seurat替换成Seurat2, seurat替换成seurat2

```bash
find . -type f -print0 | xargs -0 sed -i "s/Seurat/Seurat2/g"
find . -type f -print0 | xargs -0 sed -i "s/seurat/seurat2/g"
```

这种无差别的替换会有一个问题，会把一些这类`http://www.satijalab.org/seurat`非代码信息中的seurat替换成seurat2，不过这并不影响实际函数的使用。

将`R/seurat.R`重名为`R/seurat2.R`

```bash
mv R/seurat.R R/seurat2.R
```

之后打包修改后的文件

```bash
tar -czf Seurat2.tar.gz Seurat
```

于是我们就可以安装我们自己修改后的R包了

```r
install.packages("Seurat2.tar.gz", repos=NULL, type="source")
```

**注意**: 目前未修改测试数据集的对象，所以不能用Seurat2来运行pbmc_small的例子

由于我已经修改好了，所以你们也不需要自己再去修改了。我将其上传到GitHub，所以可以用devtools进行安装。

```r
devtools::install_github("xuzhougeng/Seurat2")
```

由于代码只是简单修改，存在bug，其实不建议尝试使用。而且安装完Seurat2之后，我**定性**的评估了下Seurat2和Seurat3，发现两者在聚类分析上并没有太大区别, 所以更推荐大家使用最新的`Seurat`.


