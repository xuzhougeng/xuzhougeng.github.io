---
title: 「小技巧」给Seurat对象瘦瘦身
date: 2022-07-04 02:33:54.935
updated: 2022-07-04 02:33:54.935
url: /archives/trick-how-to-slim-seurat-memory
categories: R
tags: 单细胞
---


我们在使用Seurat处理单细胞数据的时候，会发现Seurat对象会不断变大，一不小心就成为内存无底洞，

例如我的一个Seurat对象就占了22.3的内存空间

```r
old_size <- object.size(seu.obj)
format(old_size, units="Gb")
# 22.3 Gb
```

如果我中途需要关闭Rstudio，那么为了保证自己的工作连续性，我就需要将内存中的20多G的数据保存到磁盘上，并在下次分析加载会内存。这一来一回，考虑到磁盘的读写速度，耗时可能就需要10多分钟。

考虑到你可能要把数据上传到网盘或者GEO数据库，那么这20多G数据所需要花费的时间，就更加超出你的想象了。

那有没有办法给Seurat对象瘦身呢？ 其实很简单，因为Seurat主要是在Scale这一步，将原本的稀疏矩阵变成了普通的矩阵，同时里面的元素都是浮点型，极其占用空间。只要我们在保存数据之前先把这个归一化的矩阵给清空，就可以让Seurat一下子瘦下来。

```r
seu.obj@assays$RNA@scale.data <- matrix()
new_size <- object.size(seu.obj)
format(new_size, units="Gb")
# 5Gb
```

上面的操作，让内存占用降低到只有原来的20%左右。但是，问题来了，内存降低的代价是什么呢？代价就是，你需要对加载的数据进行scale，复原Seurat中的scale.data。

```r
all.genes <- rownames(seu.obj)
seu.obj <- ScaleData(seu.obj, features = all.genes)
```

这就是计算机科学中常见思路，要么是空间换时间，要么是时间换空间。


这个小技巧除了能加速Seurat对象的保存和读写外，还有什么其他应用吗？因为scale后数据主要是给主成分分析(PCA)提供输入。后续的非线性降维(UMAP), 聚类分析(Cluster)都是基于PCA，而不是基于scale数据，因此，如果分析过程中发现内存空间吃惊，也可以先通过这个小技巧释放下内存空间。等跑完一些占内存的操作后，再恢复即可。
