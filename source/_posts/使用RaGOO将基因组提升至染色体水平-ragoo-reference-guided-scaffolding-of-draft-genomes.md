---
title: 使用RaGOO将基因组提升至染色体水平
date: 2020-05-19 21:47:43.996
updated: 2020-05-19 21:47:43.996
url: /archives/ragoo-reference-guided-scaffolding-of-draft-genomes
categories: 生信软件工具箱 | 基因组学
tags: 组装
---

![RaGOO](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/05/image-668b73a1f1a949f1ac5b4d867cdd05d2.png)

将染色体从contig/scaffold水平提升到chromosome水平是组装的最终目标。我们通常使用遗传图谱，光学图谱，HiC这些技术提供的信息将contig进行排序连接。

如果你组装的物种刚好有一个近源物种甚至说是同一物种，那其实我们可以直接将我们的contig比对到该基因组上，根据其提供的位置信息，将我们的contig/scaffold提高到染色体水平。RaGOO就是其中一款软件，相对于其他同类型的工具，它有以下优势

- 不错的性能（感谢minimap2）
- contig完整的排序和方向调整
- GFF lift-over
- 结构变异检测，整合了Assemblytics
- 对于每个contig都计算可信得分

RaGOO使用minimap2将contig和reference进行比对，过滤低于1k的alignment，之后根据contig的覆盖度将contig聚类到最接近的染色体上，最后根据contig在染色体上的相对位置信息进行排序合并。

RaGOO基于Python3以及预先安装的minimap2。 我们需要从Github上克隆该项目进行安装

```bash
git clone https://github.com/malonge/RaGOO.git
cd RaGOO
python setup.py install
```

它的使用非常简单，就两个输入文件，contig和reference的FASTA文件

```bash
ragoo.py contigs.fasta reference.fasta
```

一些可供修改的参数:

- e: 用于忽略reference一些序列
- -gff: 将之前contig注释的GFF文件调整为当前版本
- -b: 打断chimeric contig
- -R: 提供额外的fastq/fasta序列辅助纠正错误组装
- -T: 对应-R参数提供序列的类型， sr表示short read, corr表示纠错后的long reads
- -t: 线程数
- -g: 两个contig之间的gap大小
- -s: 分析结构变异
- -i: 最低得分用于将contig分组，默认是0.2
- -j: 哪些contig序列需要忽略
- -C: 将无法锚定的contig单独成行，而非合并成一个Chr0

几个建议: 默认线程是3，可以按照自己的需求进行提高。 如果对组装没有信心，可以加上`-b -R -T`参数用来纠正潜在的错误。我强烈推荐加上`-C`, 不然你会以为`Chr0`也是一个染色体。

可能的输出文件如下

```bash
ragoo_output/
├── ctg_alignments: 错误纠正结果
├── groupings     : 分组结果
├── orderings     : 排序结果
├── pm_alns       : 结构变异分析结果
└── ragoo.fasta   : 你需要的输出文件
```

在ragoo.fasta中，默认参数下`Chr0_RaGOO`表示contig.fasta的序列无法在reference.fa中定位，直接前后相连成一个序列。

个人主观评价： RaGOO使用容易，运行效率也很高，还能够分析结构变异。根据它的文章，有些时候表现还优于HiC组装结果，以后的一些基因组项目建议用上它。

![HiC搞不定的情况](https://upload-images.jianshu.io/upload_images/2013053-50a162cda9432a11.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

[https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1829-6](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1829-6)
