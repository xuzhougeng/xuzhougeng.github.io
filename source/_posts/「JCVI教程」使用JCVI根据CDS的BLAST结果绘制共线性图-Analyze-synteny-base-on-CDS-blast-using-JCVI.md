---
title: 「JCVI教程」使用JCVI根据CDS的BLAST结果绘制共线性图
date: 2019-09-02 12:07:51.79
updated: 2019-09-02 12:07:51.79
url: /archives/Analyze-synteny-base-on-CDS-blast-using-JCVI
categories: 生信软件工具箱
tags: 可视化 | JCVI
---

> 这是唐海宝老师GitHub上的[JCVI](https://github.com/tanghaibao/jcvi)工具的非官方说明书。
> 该工具集的功能非常多，但是教程资料目前看起来并不多，因此为了能让更多人用上那么好用的工具，我就一边探索，一边写教程

这一篇文章教大家如何利用JCVI里面的工具绘制点图，展现两个物种之间的共线性关系。

在分析之前，你需要从[PhytozomeV11](https://genome.jgi.doe.gov/) 下载A.thaliana和Alyrata的CDS序列，保证文件夹里有如下内容

```bash
Alyrata_384_v2.1.cds.fa.gz        Athaliana_167_TAIR10.cds.fa.gz
Alyrata_384_v2.1.gene.gff3.gz  Athaliana_167_TAIR10.gene.gff3.gz
```

## 准备最长CDS和BED文件

我们在做CDS相互比对的时候只需要有每个基因最长的转录本即可，有两种方法可以实现

### 方法1：自己写脚本

我用我写的一个脚本`get_the_longest_transcripts.py`提取每个基因的最长转录本,见 [基因组共线性工具MCScanX使用说明](https://www.jianshu.com/p/8373e50722f6)


```bash
zcat Alyrata_384_v2.1.gene.gff3.gz | python ~/scripts/python/get_the_longest_transcripts.py > aly_lst_gene.txt
zcat Athaliana_167_TAIR10.gene.gff3.gz | python ~/scripts/python/get_the_longest_transcripts.py  > ath_lst_gene.txt
```

其中`xxx_lst_gene.txt`的格式如下, 第一列是基因名，第二列是mRNA编号，后面几列是位置信息。

```bash
$ head ath_lst_gene.txt
AT4G19470.TAIR10	AT4G19470.1.TAIR10	Chr4	10612993	10614339	-
AT5G43860.TAIR10	AT5G43860.1.TAIR10	Chr5	17630450	17632312	+
AT1G68650.TAIR10	AT1G68650.1.TAIR10	Chr1	25775741	25777874	+
AT1G28050.TAIR10	AT1G28050.1.TAIR10	Chr1	9775528	9777810	-
AT3G59880.TAIR10	AT3G59880.1.TAIR10	Chr3	22120969	22121700	+
AT1G22030.TAIR10	AT1G22030.1.TAIR10	Chr1	7759164	7760556	-
AT5G24330.TAIR10	AT5G24330.1.TAIR10	Chr5	8295147	8297068	-
AT5G43990.TAIR10	AT5G43990.2.TAIR10	Chr5	17697889	17702005	+
AT1G11410.TAIR10	AT1G11410.1.TAIR10	Chr1	3841286	3844432	+
AT4G32890.TAIR10	AT4G32890.1.TAIR10	Chr4	15875470	15876762	+
```

由于基因名和mRNA编号里有在提取CDS不需要的内容，因此要进行删除

```bash
sed -i 's/\.v2\.1//g' aly_lst_gene.txt
sed -i 's/\.TAIR10//g' ath_lst_gene.txt
```

之后我们就可以根据第二列进行提取CDS

```bash
seqkit grep -f  <(cut -f 2 ath_lst_gene.txt ) Athaliana_167_TAIR10.cds.fa.gz > ath.cds
seqkit grep -f  <(cut -f 2 aly_lst_gene.txt ) Alyrata_384_v2.1.cds.fa.gz > aly.cds
```

提取的CDS编号里面也有一些不需要的内容，所以也要删除

```bash
sed -i 's/\.t.*//' aly.cds
sed -i 's/\..*//' ath.cds
```

此外还需要基因的位置信息的bed文件

```bash
awk '{print $3"\t"$4"\t"$5"\t"$1"\t0\t"$6}' ath_lst_gene.txt | sort -k4,4V > ath.bed
awk '{print $3"\t"$4"\t"$5"\t"$1"\t0\t"$6}' aly_lst_gene.txt | sort -k4,4V > aly.bed
```

### 基于JCVI工具集

当然也可以参考[「JCVI教程」如何基于编码序列或蛋白序列进行共线性分析](https://www.jianshu.com/p/7ccc911c9273)来提取bed和cds序列，不需要用到我写的脚本。

```bash
python -m jcvi.formats.gff bed --type=mRNA --key=Name Athaliana_167_TAIR10.gene.gff3.gz  > ath.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name Alyrata_384_v2.1.gene.gff3.gz > aly.bed
```

对bed文件中的基因进行去重

```bash
python -m jcvi.formats.bed uniq ath.bed
python -m jcvi.formats.bed uniq aly.bed
```

这一步会得到`aly.uniq.bed`和`ath.uniq.bed`, 我们将其覆盖原文件

```bash
mv ath.uniq.bed ath.bed
mv aly.uniq.bed aly.bed
```

根据bed文件提取cds里的序列

```bash
# Athaliana
seqkit grep -f <(cut -f 4 ath.bed ) Athaliana_167_TAIR10.cds.fa.gz | seqkit seq -i > ath.cds
# Alyrata
seqkit grep -f <(cut -f 4 aly.bed )  Alyrata_384_v2.1.cds.fa.gz | seqkit seq -i  > aly.cds
```

---

无论是哪种方法，请保证最后有以下四个文件

```bash
$ ls ???.???
aly.bed  aly.cds  ath.bed  ath.cds
```

## BLAST比对

相对于上一步，这一步其实非常简单了

```bash
makeblastdb -in ath.cds -out db/ath -dbtype nucl
blastn -num_threads 20  -query aly.cds -db db/ath -outfmt 6 -evalue 1e-5 -num_alignments 5 > aly_ath.blast
```

用`jcvi.compara.blastfilter `对结果进行过滤

```bash
python -m jcvi.compara.blastfilter --no_strip_names aly_ath.blast --sbed ath.bed --qbed aly.bed
```

运行过程中有如下输出信息

```bash
19:59:15 [base] Load file `aly.bed`
19:59:16 [base] Load file `ath.bed`
19:59:16 [blastfilter] Load BLAST file `aly_ath.blast` (total 49887 lines)
19:59:16 [base] Load file `aly_ath.blast`
19:59:16 [blastfilter] running the cscore filter (cscore>=0.70) ..
19:59:16 [blastfilter] after filter (42023->26531) ..
19:59:16 [blastfilter] running the local dups filter (tandem_Nmax=10) ..
19:59:16 [blastfilter] after filter (26531->24242) ..
```

最后输出`aly_ath.blast.filtered`用于做图

```bash
python -m jcvi.graphics.blastplot aly_ath.blast.filtered --sbed ath.bed --qbed aly.bed
```

最后点图如下

![点图](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/9/2013053-6fbdb7b7cd866643-8edb0016df784aca86abfb4de1f7b48f.png)