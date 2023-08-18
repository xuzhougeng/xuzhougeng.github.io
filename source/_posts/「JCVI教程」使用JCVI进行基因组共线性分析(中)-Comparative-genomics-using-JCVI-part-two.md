---
title: 「JCVI教程」使用JCVI进行基因组共线性分析(中)
date: 2019-09-02 11:36:30.92
updated: 2019-09-02 11:36:30.92
url: /archives/Comparative-genomics-using-JCVI-part-two
categories: 生信软件工具箱
tags: 可视化 | JCVI
---

在[「JCVI教程」使用JCVI进行基因组共线性分析(上)](/archives/Comparative-genomics-using-JCVI-part-one)还是有一个尴尬的事情，就是只用到两个物种，不能展示出JCVI画图的方便之处，因此这里参考[https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version)](https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version))的分析，只不过画图部分拓展下思路。

首先要运行如下代码获取目的数据

```bash
python -m jcvi.apps.fetch phytozome Vvinifera,Ppersica,Tcacao
python -m jcvi.formats.gff bed --type=mRNA --key=Name Vvinifera_145_gene.gff3.gz -o grape.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name Ppersica_139_gene.gff3.gz -o peach.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name Tcacao_233_gene.gff3.gz -o cacao.bed
python -m jcvi.formats.fasta format --sep="|" Vvinifera_145_cds.fa.gz grape.cds
python -m jcvi.formats.fasta format --sep="|" Ppersica_139_cds.fa.gz peach.cds
python -m jcvi.formats.fasta format --sep="|" Tcacao_233_cds.fa.gz cacao.cds
# find ortholog
python -m jcvi.compara.catalog ortholog grape peach --cscore=.99
python -m jcvi.compara.catalog ortholog peach cacao --cscore=.99
# build .simpe
python -m jcvi.compara.synteny screen --minspan=30 --simple grape.peach.anchors grape.peach.anchors.new 
python -m jcvi.compara.synteny screen --minspan=30 --simple peach.cacao.anchors peach.cacao.anchors.new
```

然后按照教程的配置文件进行画图

layout文件内容如下

```bash
# y, xstart, xend, rotation, color, label, va,  bed
 .7,     .2,    .4,      45,      , Grape, top, grape.bed
 .5,     .2,    .8,       0,      , Peach, top, peach.bed
 .7,     .4,    .8,     -45,      , Cacao, bottom, cacao.bed
# edges
e, 0, 1, grape.peach.anchors.simple
e, 1, 2, peach.cacao.anchors.simple
```

seqids文件内容如下

```bash
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19
scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8
scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8,scaffold_9,scaffold_10r
```

运行`python -m jcvi.graphics.karyotype seqids layout`会得到如下结果

![三个物种的共线性图](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/9/2013053-9297f82212a67096-b45b77baf0e442168bef05dae6b30bb5.png)

我就在思考一个问题如何让他形成一个三角形。经过一波三角运算和不断尝试，我定义了如下的layout

```bash
# y, xstart, xend, rotation, color, label, va,  bed
 .5,      0.025,    0.625,      60,      , Grape, top, grape.bed
 .2,      0.2,    .8,       0,      , Peach, top, peach.bed
 .5,     0.375,    0.975,     -60,      , Cacao, top, cacao.bed
# edges
e, 0, 1, grape.peach.anchors.simple
e, 1, 2, peach.cacao.anchors.simple
```

那么效果怎么样呢？运行`python -m jcvi.graphics.karyotype seqids layout`吧

![三角形](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/9/2013053-3a686b648bf45bb6-602270724f8c483f8a43b19d80ae34b8.png)

当然这里只展示了，grape和peach, peach和cacao之间的共线性，我又想着能不能加上grape和caocao呢？我尝试着运行下面的代码，

```bash
python -m jcvi.compara.catalog ortholog grape cacao --cscore=.99
python -m jcvi.compara.synteny screen --minspan=30 --simple grape.cacao.anchors grape.cacao.anchors.new
```

并继续修改了layout

```bash
# y, xstart, xend, rotation, color, label, va,  bed
 .5,      0.025,    0.625,      60,      , Grape, top, grape.bed
 .2,      0.2,    .8,       0,      , Peach, top, peach.bed
 .5,     0.375,    0.975,     -60,      , Cacao, top, cacao.bed
# edges
e, 0, 1, grape.peach.anchors.simple
e, 1, 2, peach.cacao.anchors.simple
e, 0, 2, grape.cacao.anchors.simple
```

运行`python -m jcvi.graphics.karyotype seqids layout`会得到如下结果

![awesome](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/9/2013053-6577c38e1c3bd97e-b6b6f46a01a447208846c9caca202ddc.png)

我觉得给我一点时间，我也能用JCVI画出下面的图了

![amazing](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/9/2013053-2bcbafc0ba133b15-b59d4a500a284516a7d9068d502fbfa5.png)

> 图片来自于[https://science.sciencemag.org/content/345/6199/950](https://science.sciencemag.org/content/345/6199/950)
