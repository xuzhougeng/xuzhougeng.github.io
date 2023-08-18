---
title: 使用purge_haplogs处理基因组杂合区域
date: 2019-08-29 13:26:33.195
updated: 2019-12-11 10:56:49.418
url: /archives/purge_haplotigs-allelic-contig-reassignment-for-third-gen-diploid-genome-assemblies
categories: 生信软件工具箱
tags: 组装
---

FALCON和Canu的组装后会得到一个单倍型融合的基因组，用来表示二倍体基因组。之后，FALCON Unzip和Supernova这类软件进一步处理其中等位基因区域，将这部分区间进行拆分。

当基因组某些区域可能有着比较高的杂合度，这会导致基因组该区域的两个单倍型被分别组装成primary contig， 而不是一个为primary contig， 另一个是associated haplotig. 如果下游分析主要关注于单倍型，这就会导致一些问题。

那么有没有解决方案呢？其实也很好办，就是找到相似度很高的contig，将他们拆分。目前已经有一些软件可以完成类似任务，如 **HaploMerger2**, **Redundans**, 这不过这些软件主要处理二代组装结果。

 `purge_haplogs`则是最近开发，用于三代组装的基因组。它根据minimap2的比对结果，通过分析比对read的覆盖度决定谁去谁留。该工具适用于单倍型组装软件，例如 Canu, FALCON或 FALCON-Unzip primary contigs, 或者是分相后的二倍体组装(Falcon-Unzip primary contigs + haplotigs 。

它的工作流程如下图所示。一般只需要两个输入文件，组装草图(FASTA格式) 和 比对的BAM文件。同时还可以提供重复序列注释的BED文件，辅助处理高重复的contig。

![分析流程](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/2013053-aa0d1c1178c3046a-0cd9a3bad7ce4e4f9c775ba526caaf67.png)

> **建议**: 用原来用于组装的read进行比对。对于多个匹配的read，建议采取random best，也就是随便找一个。

## 软件安装

`purge_haplotigs`依赖软件比较多，手动安装会很麻烦，但是他可以直接用bioconda装

```bash
conda create -n purge_haplotigs_env
conda activate purge_haplotigs_env
conda install purge_haplotigs
```

安装完成后需要一步测试

```bash
purge_haplotigs test
```

## 简明教程

数据准备。 需要下载的数据集分为两个部分，一个是FALCON-Unzip后的primary contig 和 halplotigs. 另一个则是已经比完后的BAM文件

```bash
mkdir purge_haplotigs_tutorial
cd purge_haplotigs_tutorial
wget https://zenodo.org/record/841398/files/cns_h_ctg.fasta
wget https://zenodo.org/record/841398/files/cns_p_ctg.aligned.sd.bam # 1.7G
wget https://zenodo.org/record/841398/files/cns_p_ctg.aligned.sd.bam.bai 
 wget https://zenodo.org/record/841398/files/cns_p_ctg.fasta
wget https://zenodo.org/record/841398/files/cns_p_ctg.fasta.fai
```

当然我们不可能直接就拿到比对好的BAM文件，我们一般是有组装后的基因组以及用于组装的subread，假设这两个文件命名为, genome.fa 和 subreads.fasta.gz.

官方提供的新比对代码

```bash
minimap2 -t 4 -ax map-pb genome.fa subreads.fasta.gz --secondary=no \
    | samtools sort -m 1G -o aligned.bam -T tmp.ali
```

如下是旧版代码

```bash
minimap2 -ax map-pb genome.fa subreads.fasta.gz \
    | samtools view -hF 256 - \
    | samtools sort -@ 8 -m 1G -o aligned.bam -T tmp.ali
```

如果你有二代测序数据，也可以用BWA-MEM进行比对得到BAM文件。

第一步：使用`purge_haplotigs  readhist`从BAM中统计read深度，绘制柱状图。

```bash
# 新
purge_haplotigs  hist  -b aligned.bam  -g genome.fasta  -t 20
# 旧
# purge_haplotigs  readhist  -b aligned.bam  -g genome.fasta  -t 20
# -t 线程数, 不宜过高，过高反而没有效果。
```

也就是下图，你能明显的看到图中有两个峰，一个是单倍型的覆盖度，另一个二倍型的覆盖度，

![高杂合基因组read-depth histogram](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/2013053-b3c0c46eba030da3-f23b4f0d95094a4eb78d42c7e6f8a875.png)

你可能还想知道高纯合基因组是什么样的效果，我也找了一个纯合的物种做了也做了read-depth 柱状图，

![纯合基因组read-depth histogram](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/2013053-77166f97f8b14c91-4965aa2fbcaf48a590a98091d9f7d0af.png)


之后你需要根据read-depth 柱状图 确定这两个峰的位置用于下一步。下面是两个例子。对于我们则是，20，65，190.

![两个例子](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/2013053-7561308e761d13ee-12e97f73e8c84bcb8e3b7051a08ceda7.png)

第二步: 根据read-depth信息选择阈值。

```bash
# 新
purge_haplotigs  contigcov  -i cns_p_ctg.aligned.sd.bam.gencov  -o coverage_stats.csv  -l 20  -m 75  -h 190
# 旧
# purge_haplotigs  contigcov  -i cns_p_ctg.aligned.sd.bam.gencov  -o coverage_stats.csv  -l 20  -m 75  -h 190
```

这一步生成的文件是"coverage_stats.csv"

第三步：区分haplotigs. 

```bash
purge_haplotigs purge  -g cns_p_ctg.fasta  -c coverage_stats.csv  -b cns_p_ctg.aligned.sd.bam  -t 4  -a 60
```

这一步会得到如下文件

- curated.artefacts.fasta：无用的contig，也就是没有足够覆盖度的contig.
- curated.fasta：新的单倍型组装
- curated.haplotigs.fasta：从原本组装分出来的haplotigs
- curated.reassignments.tsv: 单倍型的分配信息
- curated.contig_associations.log: 运行日志, 下面是其中一个记录，表示000004F_004和000004F_027是000004F_017的HAPLOTIG, 而000004F_017和000004F_013又是000004F,的HAPLOTIG。

```bash
000004F,PRIMARY -> 000004F_013,HAPLOTIG
                -> 000004F_017,HAPLOTIG 
                                        -> 000004F_004,HAPLOTIG
                                        -> 000004F_027,HAPLOTIG
```

在"curated.reassignments.tsv"文件中有6列

- reassigned_contig: 用于比较的contig
- top_hit_contig: 最好的被比对的contig
- second_hit_contig: 第二个被比对的contig
- best_match_coverage: 最好的匹配覆盖度
- max_match_coverage :  最高的匹配深度
- reassignment: 标记为haplotype 还是 repeat，或者是keep

由于我们用的是单倍型组装primary contigs而不是二倍体组装的parimary + haplotigs, 因此我们需要将FALCON_Unzip的haplotgi合并到重新分配的haplotigs中，这样子我们依旧拥有二倍体组装结果

```bash
cat cns_h_ctg.fasta >> curated.haplotigs.fasta
```

### 检查dotplots

如果在第三步`purge_haplotigs purge`中添加了`-d/--dotplots`参数，即为每个reassigned_contigs和unassigned_contigs生成用于人工检查的共线性图，那么在最终的输出结果中会有两个文件目录

- dotplots_reassigned_contigs
- dotplots_unassigned_contigs

官方提供了5个可能出现的共线性图

第一种: [Haplotig](https://bitbucket.org/repo/Ej8Mz7/images/2663721414-haplotig.png)，最佳情况，完美的共线性关系

![Haplotig](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/12/image-20191211103643627-2242c4d3dbb646019a52ec5d23c3a4a7.png)

第二种: [大部分是haplotigs](https://bitbucket.org/repo/Ej8Mz7/images/3628840007-haploid_diploid_hemizygous.png), 说明这个contig部分是二倍型，部分是单倍型，可能是半合子(hemizygous)

![mostly haplotig](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/12/image-20191211103820925-4dc398c1eef94efa8e6fe42c49ae8dd0.png)

第三种: Haplotigs里有大量的[串联重复](https://bitbucket.org/repo/Ej8Mz7/images/854324485-repeat_rich.png)

![tandem repeat](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/12/image-20191211104431607-4df146fcbf1d4e178c13c2d205430204.png)

第四种: Haplotigs是[回文序列(palindrome)](https://bitbucket.org/repo/Ej8Mz7/images/1079545124-haplotig_with_palindrome.png)

![palindrome](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/12/image-20191211104300406-d4ff89a584c64324812cb1c3dd407d56.png)

第五种: contig从string graph中[knots](https://bitbucket.org/repo/Ej8Mz7/images/228451888-repeats_string_graph_short_cut.png)产生，这种情况不算是haplotigs，但是对于短读序列的比对会造成麻烦

![knots](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/12/image-20191211104509706-8895b6b368bc4255aa702af868f3e7af.png)

你可能需要看大量的图才能有感觉，到底应该把哪些Purge_haplotigs错误认为是haplotig的contig放回到primary contig中。


## 原理介绍

为什么第一步的 read 覆盖深度分析能判断基因组是否冗余呢？这是因为对于坍缩的单倍型，那么含有等位的基因的read只能比对到该位置上，而如果杂合度太高被拆分成两个不同的contig，那么含有等位的基因的read就会分别比对到不同的read上，导致深度降低一半。下图A就是一个典型的包含冗余基因组的read覆盖度分布

![read-深度分析](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/2013053-2f26047d1da456f5-72ba2aa4116245759f400928dae6d383.png)


分析流程的第二步的任务就是人工划分出如下图B部分，绿色的部分是坍缩单倍型contig，蓝色的部分是潜在的冗余contig。之后，分析流程会计算这些区域中的contig的覆盖度。 对于绿色部分中的contig，如果覆盖度低于80%, 会进行标记用于后续分析。如果深度非常低，那么很有可能就是组装引入错误，深度非常高的部分基本就是重复序列或者是细胞器的contig，这些黄色的contig可以在后续的组装出分开。

![划分区间](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/2013053-0ae4a373254ca9e8-a7ccd5ff06a94fd08c5001b57f6c1f52.png)

第三步就是同源序列进行识别和分配。所有标记的contig之后会用Minimap2在整个组装进行搜索，寻找相似度较高的离散区间（如下图C）。如果一个Contig的联配得分大于阈值(默认70%), 那么就会被标记为haplotigs. 如果一个contig的最大联配得分大于阈值(默认250%), 会被标记成重复序列，这有可能是潜在的有问题contig，或许是坍缩的contig或者低复杂度序列。


![移除haplotigs](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/2013053-ca64a560ff4099bd-c42225d428ee40bb9f5feb048e7963a9.png)


## 推荐阅读

-  Purge Haplotigs: allelic contig reassignment for third-gen diploid genome assemblies
-  <https://bitbucket.org/mroachawri/purge_haplotigs>