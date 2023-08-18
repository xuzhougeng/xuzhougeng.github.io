---
title: 如何使用GMAP/GSNAP进行转录组序列比对
date: 2022-01-18 08:06:16.614
updated: 2022-01-18 08:06:16.614
url: /archives/align-transcript-with-gmap
categories: 转录组学
tags: 注释
---

GMAP最早用于讲EST/cDNA序列比对到参考基因组上，可以用于基因组结构注释。后来高通量测序时代，又开发了GSNAP支持高通量数据比对，这篇文章主要介绍GMAP，毕竟高通量转录组数据比对大家更喜欢用STAR, HISTA2等软件。

## 软件安装

下面是我源码安装的代码

```bash
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2018-07-04.tar.gz
tar xf gmap-gsnap-2018-07-04.tar.gz
cd gmap-2018-07-04/
./configure --prefix=$HOME/opt/biosoft/gmap
make -j 20
```

## 软件使用

如下步骤假设你有一个物种的基因组序列和对应的CDS序列，分别命名为"reference.fa"和"cds.fa"

### 第一步:构建GMAP/GSNAP索引数据库

GMAP/GSNAP对FASTA文件中每个记录下的序列的长度有一定限制, 每一条不能超过4G, 能应付的了大部分物种了。

构建索引分为两种情况考虑，第一种是一个fasta文件包含所有的序列

```bash
~/opt/biosoft/gmap/bin/gmap_build -d reference reference.fa
```

第二种则是每个染色体的序列都单独存放在一个文件夹里，比如说你下载人类参考基因组序列解压后发现有N多个fasta文件, 然后你就想用其中几条染色体构建索引

```bash
~/opt/biosoft/gmap/bin/gmap_build -d reference Chr1.fa Chr2.fa Chr3.fa ...
```

注: 这里的`-d`表示数据K库的名字，默认把索引存放在gmap安装路径下的share里，可以用`-D`更改.此外还有一个参数`-k`用于设置K-mer的长度, 默认是15, 理论上只有大于4GB基因组才会有两条一摸一样的15bp序列(当然是完全随机情况下)。

### 第二步:正式比对

建立完索引之后就可以将已有的CDS或者EST序列和参考基因组序列进行比较。

```bash
~/opt/biosoft/gmap/bin/gmap -t 10 -d reference -f gff3_gene cds.fa > cds_gene.gff3
```

其中`-t`设置线程数, `-d`表示参考基因组数据库的名字, 都是常规参数。我比较感兴趣的参数是如何将序列输出成GFF格式. GMAP允许多种格式的输出，比如说`-S`只看联配的总体情况，而`-A`会显示每个比对上序列的联配情况, 还可以输出蛋白序列(-P)或者是genomic序列(-E). 但是做结构注释要的gff文件，参数就是`-f gff3_gene`, `-f gff3_match_cdna`, `-f gff3_match_est`。

## 参考文献

要想对一个软件有更好的认识，最好还是看看他们文章是怎么说的。

- GMAP: a genomic mapping and alignment program for mRNA and EST sequences 
Bioinformatics 2005 21:1859-1875 [Abstract](http://bioinformatics.oupjournals.org/cgi/content/abstract/21/9/1859) [Full Text](http://bioinformatics.oupjournals.org/cgi/content/full/21/9/1859), Thomas D. Wu and Colin K. Watanabe
- Fast and SNP-tolerant detection of complex variants and splicing in short reads 
Bioinformatics 2010 26:873-881 [Abstract](http://bioinformatics.oupjournals.org/cgi/content/abstract/26/7/873)[Full Text](http://bioinformatics.oupjournals.org/cgi/content/full/26/7/873), Thomas D. Wu and Serban Nacu
 
