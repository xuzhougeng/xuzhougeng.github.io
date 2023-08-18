---
title: 三代组装软件miniasm笔记
date: 2019-08-29 22:50:12.814
updated: 2019-08-29 22:51:25.896
url: /archives/Assemble-nanopore-using-minimap-and-miniasm
categories: 生信软件工具箱
tags: 组装
---

我们用来练手的文章发表在 _Nature Communication_ ，"High contiguity Arabidopsis thaliana genome assembly with a single nanopore flow cell", 非常不要脸的说，这篇文章是我师爷实验室发的。

简单讲讲故事内容，就是他们实验室买了一台nanopore仪器，就是下面这台， 目前仪器价格国内是8K左右，当然测序的价格就另说了。如同买台PS4主机，还要买游戏，买个单反，你还得买镜头。仪器只是败家的开始！

![nanopore](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/2013053-ada4497b7b113f99-90a7dabd172345d7bfb5087f9b4b9587.png)

他们认为三代测序目前有两大问题，测的还不够长以及不够准。nanopore解决了其中一个问题，不够长。_Arabidopsis thaliana_ 当年用一代测序，虽然可以认为是组装的金标准了，但是还是有很多区域是BAC连BAC文库搞不定的，所以就用这台仪器把 _Arabidopsis thaliana_ 测了一波。显然就测一个nanopore，还是已知序列的物种是不可能发文章的，于是他们又用Pacbio sequel测了一波。最后用bionano 光学图谱验证了一次(请大家自行计算要多少钱)。

光测序不行，还得组装对吧。传统的组装方法是想办法利用高深度和随机错误进行纠错，然后用纠错后的长序列进行组装，最后用二代进行纠错。对于一台不错的服务器（20W起步吧）大约花个十天半个月就行。作者或许认为买一台20多w的外设配合不到1w的测序仪可能是太蠢了，于是他用了比较Li Heng大神开发的工具，Minimap+miniasm进行组装，然后用racon+pillon进行纠错，用了一台Macbook Pro 15.6寸花了4天就搞定了，并且和常规工具比较，还算过得去哦。

下面就是正式的分析：

根据文章提供的项目编号"PRJEB21270", 在European Nucleotide Archive上找到下载地址。

![ENA搜索](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/2013053-57cf8b9704a5c90d-154ad858607d49858f55d5ab0905be20.png)

进入这个页面之后，就可以去下载作者用到的所有数据，我们下载Sequel和MinIon和Illuminia的数据就好了，数据量加起来差不多30G。

```bash
## Sequal
wget -c -q ftp://ftp.sra.ebi.ac.uk/vol1/ERA111/ERA1116568/bam/pb.bam
wget -c -q ftp://ftp.sra.ebi.ac.uk/vol1/ERA111/ERA1116568/bam/pb.bam.bai
## MinION
wget -c -q ftp://ftp.sra.ebi.ac.uk/vol1/ERA111/ERA1116595/fastq/ont.fq.gz
```

对于Illumina的二代测序，需要用prefetch进行下载

```bash
# Illuminia MiSeq
prefetch ERR2173372
fasterq-dump -O  . ERR2173372
```

拿到数据之后，我们就可以用作者提供的分析流程进行重复了。地址为<https://github.com/fbemm/onefc-oneasm/wiki/Assembly-Generation>

> 这就是大神的自信，把代码都给你，反正你也看不懂。当然我在重复的时候用的都是最新的软件，所以会有所不同

第一步：拿着80%～90%正确率的原始数据相互比对， 找序列之间的Overlap。这一步，我花了30分钟

```bash
time ~/opt/biosoft/minimap2/minimap2 -t 10 -x ava-ont ont.fq ont.fq > gzip -1 ont.paf.gz &
```

第二步：找到Overlap，就能够进行组装了。这一步我花了2分钟

```bash
time ~/opt/biosoft/miniasm/miniasm -f ont.fq ont.paf > ONTmin.gfa &
awk '/^S/{print ">"$2"\n"$3}' ONTmin.gfa | seqkit seq > ONTmin_IT0.fasta &
```

第三步： 原始的组装结果充满了错误，所以需要进行纠错。纠错分为两种，一种是用三代自身数据，一种是用二代数据进行纠错。当然这两步都是需要的

首先使用三代数据进行纠错，古语有云“事不过三”一般迭代个三次就差不多。这三步，差不多用了1个小时。

```bash
# Iteration 1
~/opt/biosoft/minimap2/minimap2 ONTmin_IT0.fasta ont.fq > ONTmin_IT0.paf &
time ~/opt/biosoft/racon/build/bin/racon -t 10 ont.fq ONTmin_IT0.paf ONTmin_IT0.fasta > ONTmin_IT1.fasta &
# Iteration 2
~/opt/biosoft/minimap2/minimap2 ONTmin_IT1.fasta ont.fq > ONTmin_IT1.paf
time ~/opt/biosoft/racon/build/bin/racon -t 10 ont.fq ONTmin_IT1.paf ONTmin_IT1.fasta> ONTmin_IT2.fasta
# Iteration 3
~/opt/biosoft/minimap2/minimap2 ONTmin_IT2.fasta ont.fq > ONTmin_IT2.paf
time ~/opt/biosoft/racon/build/bin/racon -t 10 ont.fq ONTmin_IT2.paf ONTmin_IT2.fasta > ONTmin_IT3.fasta
```

之后使用二代数据进行纠错。二代数据虽然短，但是测序质量高，所以一般都要用它进行纠错。推荐用30X PCR free的illuminia 测序数据。

Step 1: 数据预处理，过滤低质量短读，去接头。工具很多，常用的是trimmomatic，cutadapter.  我安利一个国内海普洛斯搞的一个工具fastp。

```bash
# data clean
fastp -q 30 -5 -l 100 -i ERR2173372_1.fastq -I ERR2173372_2.fastq -o i1_clean_1.fq -O i1_clean_2.fq 
```

这里标准为：平均质量高于Q30，对5‘端进行低质量碱基删除，保留大于100bp的短读

Step2:  比对，这一步基本都只用了bwa了

```bash
# align
bwa index ONTmin_IT3.fasta
bwa mem -t 8 ONTmin_IT3.fasta il_clean_1.fastq il_clean_2.fastq | samtools sort -@ 8 > ONTmin_IT3.bam
```

step3: 使用比对后的BAM文件进行纠错

```bash
# short read consensus call
java -Xmx16G -jar pilon-1.22.jar --genome ONTmin_IT3.fasta --frags ONTmin_IT3.bam --fix snps --output ONTmin_IT4
```

二代纠错的时间明显比之前的久，需要一天时间。

大家拿出自己的笔记本实际感受下呗

## 参考文献

- nanopore组装拟南芥: High contiguity Arabidopsis thaliana genome assembly with a single nanopore flow cell
- 不纠错组装: Minimap and miniasm: fast mapping and de novo assembly for noisy long sequences
- 三代组装软件评测: Comprehensive evaluation of non-hybrid genome assembly tools for third-generation PacBio long-read sequence data