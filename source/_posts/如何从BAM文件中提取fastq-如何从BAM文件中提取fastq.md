---
title: 如何从BAM文件中提取fastq
date: 2019-07-29 10:03:14.149
updated: 2019-07-29 10:04:30.371
url: /archives/如何从BAM文件中提取fastq
categories: Linux
tags: 
---

虽然高通量测序分析最常用的操作是将fastq比对到参考基因组得到BAM文件，但偶尔我们也需要提取BAM文件中特定区域中fastq。最开始我认为这是一个非常简单的操作，因为samtools其实已经提供了相应的工具`samtools fastq`. 

以biostar handbook的Ebola病毒数据为例，首先获取比对得到的BAM文件。

```bash
# 建立文件夹
mkdir -p refs
# 根据Accession下载
ACC=AF086833
REF=refs/$ACC.fa
efetch -db=nuccore -format=fasta -id=$ACC | seqret -filter -sid $ACC > $REF
bwa index $REF
SRR=SRR1553500
fastq-dump -X 100000 --split-files $SRR
BAM=$SRR.bam
R1=${SRR}_1.fastq
R2=${SRR}_2.fastq
TAG="@RG\tID:$SRR\tSM:$SRR\tLB:$SRR"
bwa mem -R $TAG $REF $R1 $R2 | samtools sort > $BAM
samtools index $BAM
```

让我们尝试提取下前1kb比对的fastq文件

```bash
samtools view -h SRR1553500.bam KJ660346.1:1-1000 | samtools fastq -1 read_1.fq -2 read_2.fq -
#[M::bam2fq_mainloop] discarded 0 singletons
#[M::bam2fq_mainloop] processed 8702 reads
```

看起来似乎我们成功了，那让我们把这些序列比对回去吧

```
bwa mem refs/AF086833.fa read_1.fq read_2.fq
[mem_sam_pe] paired reads have different names: "SRR1553500.11772", "SRR1553500.43775"
```

你会发现这里产生了错误。这个错误说的是配对的read名字不同。导致问题的原因是：samtools fastq会按照fastq在bam里原本的顺序进行提取，而按照位置排序的bam的两个配对read并不是前后紧挨着的。所以我们需要进行一波排序。

```bash
cat read_1.fq \
| paste - - - - \
| sort -k1,1 -S 3G \
| tr '\t' '\n' \
| gzip > read_1.fq.gz
cat read_2.fq \
| paste - - - - \
| sort -k1,1 -S 3G \
| tr '\t' '\n' \
| gzip > read_2.fq.gz
```

这里将5个shell命令用管道连起来完成这个任务，咋看起来比较复杂，但是理解每一步的目的，每个参数的作用，以后就能活用起来了。

首先用cat**逐行**读取read，传入到paste中。paste的作用是将多个文件合并成一个文件，做法比较简单粗暴，举个例子

```bash
cat A.txt
A
B
C
cat B.txt
D
E
F
paste A.txt B.txt
A	D
B	E
C	F
```

`cat read_1.fq | paste - - - - `表示读取4行合并成一行。相当于fastq里的read都会用一行而不是四行表示。接着用`sort`按照第一行(-k1,1)排序，内存的缓冲数据为3G(-S 3G). 然后将制表符替换成换行符，让一行又变回到4行，后续是压缩成gz，降低磁盘占用。

这样做其实还没又解决问题，因为这两个文件里的reads数其实不相同，我们需要进一步从里面剔除掉只有一方才有的read。这里使用的`seqkit`，一个非常强大的序列处理轮子。

```bash
seqkit grep -f <(seqkit seq -n -i read_1.fq.gz ) read_2.fq.gz  > read_flt_2.fq
seqkit grep -f <(seqkit seq -n -i read_2.fq.gz ) read_1.fq.gz  > read_flt_1.fq
wc -l read_flt_1.fq
   15028 read_flt_1.fq
wc -l read_flt_2.fq
   15028 read_flt_2.fq
```

终于，你得到了配对存放的两个fastq文件。

通过这个案例我想表达的观点是：虽然有些时候没有现成的工具可以解决你的问题，但是你可以通过组合几个现有的工具来完成。也就是将一个大问题拆解成几个独立的小问题，然后逐个击破。

当然，你首先得熟悉一些已有的工具，比如说常用unix命令:sort, paste, tr等，和一些非常强大生信命令行工具：samtools, bedtools, seqkit等。

## 一个更佳方便的方法

前面只是问题的一种方法，可能还不是最好的，但是我当初最先想到的一个比较粗暴的方法。当然更佳巧妙的方法是，提取序列之后可以按照read name排序，然后提取

```bash
samtools view -h SRR1553500.bam KJ660346.1:1-1000 | samtools sort -n | samtools fastq -1 read_1.fq -2 read_2.fq -s singleton.fq -
```

是不是更好理解了。


