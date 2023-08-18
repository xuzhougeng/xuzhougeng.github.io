---
title: 「热图」ComplexHeatmap展示单细胞聚类
date: 2019-10-24 19:49:27.742
updated: 2019-10-24 21:15:12.426
url: /archives/ComplexHeatmap-for-single-cell-visualization
categories: R
tags: 可视化 | 单细胞
---


实用Seurat自带的热图函数`DoHeatmap`绘制的热图，感觉有点不上档次，于是我尝试使用`ComplexHeatmap`这个R包来对结果进行展示。

个人觉得好的热图有三个要素

- 聚类: 能够让别人一眼就看到模式
- 注释: 附加注释能提供更多信息
- 配色: 要符合直觉，比如说大部分都会认为红色是高表达，蓝色是低表达

在正式开始之前，我们需要先获取一下pbmc的数据，Seurat提供了R包`SeuratData`专门用于获取数据

```r
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
InstallData("pbmc3k")
```

加载数据并进行数据预处理，获取绘制热图所需的数据

```r
library(SeuratData)
library(Seurat)
data("pbmc3k")
pbmc <- pbmc3k
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)
```

先感受下Seurat自带热图

```r
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
```

![Seurat-heatmap](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/10/1571906906510-5666cad1b38f473b8b8896f9ebde28d3.png)

下面则是介绍如何用R包`ComplexHeatmap`进行组图，虽然这个R包名带着Complex，但是并不是说这个R包很复杂，这个Complex应该翻译成复合，也就是说这个R包能在热图的基础上整合很多信息。

先安装并加载R包。

```r
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
```

为了手动绘制一个热图，要从Seurat对象中提取所需要的表达量矩阵。我提取的是原始的count值，然后用`log2(count + 1)`的方式进行标准化

```R
mat <- GetAssayData(pbmc, slot = "counts")
mat <- log2(mat + 1)
```

获取基因和细胞聚类信息

```r
gene_features <- top10
cluster_info <- sort(pbmc$seurat_annotations)
```

对表达量矩阵进行排序和筛选

```r
mat <- as.matrix(mat[top10$gene, names(cluster_info)])
```

用`Heatmap`绘制热图。对于单细胞这种数据，一定要设置如下4个参数

- `cluster_rows= FALSE`: 不作行聚类
- `cluster_columns= FALSE`: 不作列聚类
- `show_column_names=FALSE`: 不展示列名
- `show_row_names=FALSE`: 不展示行名，基因数目不多时候可以考虑设置为TRUE

```r
Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE)
```

![Heatmap-1](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/10/1571911526690-f312b4bf617f4e68bb539556872fe405.png)

从图中，我们可以发现以下几个问题：

- 长宽比不合理，当然这和绘图函数无关，可以在保存时修改长宽比
- 基因名重叠，考虑调整大小，或者不展示，或者只展示重要的基因
- 颜色可以调整
- 缺少聚类信息

这些问题，我们可以通过在[ComplexHeatmap Complete Reference](https://jokergoo.github.io/ComplexHeatmap-reference/book/)查找对应信息来解决。

## 配色方案

在热图中会涉及到两类配色，一种用来表示表达量的连续性变化，一种则是展示聚类。有一个神奇的R包就是用于处理配色，他的Github地址为<https://github.com/caleblareau/BuenColors>。

```bash
devtools::install_github("caleblareau/BuenColors")
library("BuenColors")
```

它提供了一些列预设的颜色，比方说`jdb_color_maps`

```r
      HSC       MPP      LMPP       CMP       CLP       MEP       GMP 
"#00441B" "#46A040" "#00AF99" "#FFC179" "#98D9E9" "#F6313E" "#FFA300" 
      pDC      mono     GMP-A     GMP-B     GMP-C       Ery       CD4 
"#C390D4" "#FF5A00" "#AFAFAF" "#7D7D7D" "#4B4B4B" "#8F1336" "#0081C9" 
      CD8        NK         B 
"#001588" "#490C65" "#BA7FD0"
```

这些颜色就能用于命名单细胞的类群，比如说我选择了前9个

```r
col <- jdb_color_maps[1:9]
names(col) <- levels(cluster_info)
```

## 增加列聚类信息

`Heatmap`的`row_split`和`column_split`参数可以通过设置分类变量对热图进行分隔。更多对热图进行拆分，可以参考[Heatmap split](https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#heatmap-split)

```r
Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info)
```

![Heatmap-2](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/10/1571911799187-41f8f8f51c7f4ee297fe3e29238daad1.png)

只用文字描述可能不够好看，最好是带有颜色的分块图，其中里面的颜色和t-SNE或UMAP聚类颜色一致，才能更好的展示信息。

为了增加聚类注释，我们需要用到`HeatmapAnnotation`函数，它对细胞的列进行注释，而`rowAnnotation`函数可以对行进行注释。这两个函数能够增加各种类型的注释，包括条形图，点图，折线图，箱线图，密度图等等，这些函数的特征是`anno_xxx`，例如`anno_block`就用来绘制区块图。

```r
top_anno <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = col), # 设置填充色
                       labels = levels(cluster_info), 
                       labels_gp = gpar(cex = 0.5, col = "white"))) # 设置字体
```

其中`anno_block`中的`gp`参数用于设置各类**图形参数**，`labels`设置标签，`labels_gp`设置和标签相关的**图形参数**。可以用`?gp`来了解有哪些**图形参数**。

```r
Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info,
        top_annotation = top_anno, # 在热图上边增加注释
        column_title = NULL ) # 不需要列标题
```

![Heatmap-3](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/10/1571914037258-69b71edb4e354f759b9f07674e714f13.png)

## 突出重要基因

由于基因很多直接展示出来，根本看不清，我们可以强调几个标记基因。用到两个函数是`rowAnnotation`和`anno_mark`

已知不同类群的标记基因如下

| Cluster ID | Markers       | Cell Type    |
| :--------- | :------------ | :----------- |
| 0          | IL7R, CCR7    | Naive CD4+ T |
| 1          | IL7R, S100A4  | Memory CD4+  |
| 2          | CD14, LYZ     | CD14+ Mono   |
| 3          | MS4A1         | B            |
| 4          | CD8A          | CD8+ T       |
| 5          | FCGR3A, MS4A7 | FCGR3A+ Mono |
| 6          | GNLY, NKG7    | NK           |
| 7          | FCER1A, CST3  | DC           |
| 8          | PPBP          | Platelet     |

我们需要给`anno_mark`提供基因所在行即可。

```r
mark_gene <- c("IL7R","CCR7","IL7R","S100A4","CD14","LYZ","MS4A1","CD8A","FCGR3A","MS4A7","GNLY","NKG7","FCER1A", "CST3","PPBP")
gene_pos <- which(rownames(mat) %in% mark_gene)

row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, 
                                     labels = mark_gene))

```

接着绘制热图

```r
Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info,
        top_annotation = top_anno,
        right_annotation = row_anno,
        column_title = NULL)
```

![Heatmap-4](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/10/1571914905222-c332ebab06ee4100838180f2cb0f698a.png)

关于如何增加标记注释，参考[mark-annotation](https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#mark-annotation)

## 调增图例位置

目前的热图还有一个问题，也就是表示表达量范围的图例太占位置了，有两种解决方法

- 参数设置`show_heatmap_legend=FALSE`直接删掉。
- 利用`heatmap_legend_param`参数更改样式

我们根据[legends](https://jokergoo.github.io/ComplexHeatmap-reference/book/legends.html)这一节的内容进行一些调整

```r
Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info,
        top_annotation = top_anno,
        right_annotation = row_anno,
        column_title = NULL,
        heatmap_legend_param = list(
          title = "log2(count+1)",
          title_position = "leftcenter-rot"
        ))
```

![heatmap-5](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/10/1571915704401-c349679da08d48918a452ba49adbd6cb.png)

因为ComplextHeatmap是基于Grid图形系统，因此可以先绘制热图，然后再用`grid::draw`绘制图例，从而实现将条形图的位置移动到图中的任意位置。

先获取绘制热图的对象

```r
p <- Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info,
        top_annotation = top_anno,
        right_annotation = row_anno,
        column_title = NULL,
        show_heatmap_legend = FALSE
        )
```

根据`p@matrix_color_mapping`获取图例的颜色的设置，然后用`Legend`构建图例

```r
col_fun  <- circlize::colorRamp2(c(0, 1, 2 ,3, 4),
                                 c("#0000FFFF", "#9A70FBFF", "#D8C6F3FF", "#FFC8B9FF", "#FF7D5DFF"))
lgd <-  Legend(col_fun = col_fun, 
               title = "log2(count+1)", 
               title_gp = gpar(col="white", cex = 0.75),
               title_position = "leftcenter-rot",
               #direction = "horizontal"
               at = c(0, 1, 4), 
               labels = c("low", "median", "high"),
               labels_gp = gpar(col="white")
               )
```

绘制图形

```r
grid.newpage() #新建画布
draw(p) # 绘制热图
draw(lgd, x = unit(0.05, "npc"), 
     y = unit(0.05, "npc"), 
     just = c("left", "bottom")) # 绘制图形
```

![heatmap-6](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/10/1571916860482-418dc81dde834ac3a2ab4e65fc65a005.png)

ComplexHeatmap绘制热图非常强大的工具，大部分我想要的功能它都有，甚至我没有想到的它也有，这个教程只是展示其中一小部分功能而已，还有很多功能要慢慢探索。

