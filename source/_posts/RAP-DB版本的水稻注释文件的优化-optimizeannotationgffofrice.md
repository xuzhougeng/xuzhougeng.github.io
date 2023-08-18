---
title: RAP-DB版本的水稻注释文件的优化
date: 2022-02-14 07:17:21.226
updated: 2022-02-17 01:52:35.556
url: /archives/optimizeannotationgffofrice
categories: 基因组学
tags: 小技巧
---

> 代码和数据都在GitHub上，见 <https://github.com/xuzhougeng/rice_annotation>

水稻有两个主要的注释项目组

- [RAP-DB](https://rapdb.dna.affrc.go.jp): 编号格式形如Os01g0100100
- [RGAP](http://rice.uga.edu):  编号格式形如LOC_Os0101010

虽然这两个注释项目组使用的是同一个水稻基因组，但是注释结果却存在一些差异。你可以使用任意一套注释进行分析，通过ID的对应关系，将其中一套编号转换成另外一套编号，我用的是RAP-DB版本，这是因为KEGG数据库使用的就是RAP-DB的编号形式。

RAP-DB的网站上提供了很多数据方便我们使用，比如说光基因序列文件就至少提供了为UTR+CDS, UTR+CDS+Intron, CDS三类，非常贴心。美中不足的一点是，他们提供的GFF文件并不规整。

他们提供了2个下载链接，分别对应GFF和GTF

- [Gene structure and function information in GFF format](https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_2021-11-11.tar.gz).
- [Gene structure (only exon) information in GTF format](https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_transcript_exon_2021-11-11.gtf.gz).


GFF下载之后是一个压缩包，解压能得到locus.gff, transcripts_exon.gff, transcripts.gff. GTF下载之后只有一个gtf文件。

然而，无论是GFF还是GTF文件，他们都没有提供完整信息。GTF里面只有transcript和exon. GFF要么只有gene(locus.gff), 要么只有mRNA和exon(transcripts_exon.gff), 要么是mRNA,UTR和CDS(transcripts.gff). 我就非常纳闷了，为什么没有一个文件，同时有gene, mRNA, UTR, CDS, exon这些完整信息呢？

算了，既然他们没有给，我就自己动手，于是我写了一些代码，基于他们提供loucs.gff和transcripts.gff进行整理，代码我上传到GitHub。

我在代码中，将chr01-12改成了Chr1-Chr12, 所以对应的基因组文件，我也手动调整了。整理好的FASTA, GTF, GFF文件，我都一并上传到GitHub上，方便后续使用，地址为 <https://github.com/xuzhougeng/rice_annotation>

既然整理好了FASTA和GFF文件，为了方便后续在R语言里面使用，我还生成了对应BSgenome和TxDb，安装方法如下（生成代码和输出R包都在GitHub上）

安装方式(以BSgenome为例):

1. 在R里安装：`install.pakcages("/路径/到/BSgenome.Osativa.MSU.xzg_1.0.tar.gz", repos=NULL, type="source) `
1. 在命令行里安装: `R CMD INSTALL BSgenome.Osativa.MSU.xzg_1.0.tar.gz`

为了避免和Biocondutor其他MSU的包冲突，所以我用自己名字的缩写xzg作为后缀。

安装之后，就可以通过library进行加载

```r
library(BSgenome.Osativa.MSU.xzg)
library(TxDb.Osativa.MSU.xzg)
```


