---
title: ALLHiC续:如何构建Allele.ctg.table
date: 2019-10-30 23:25:18.251
updated: 2019-10-30 23:48:27.104
url: /archives/How-to-build-allele-ctg-table
categories: 生信软件工具箱
tags: 组装 | Hi-C
---


Allele.ctg.table是[ALLHiC]( /archives/Assemble-genome-using-ALLHiC-with-HiC-Data )用于过滤多倍体基因组中因等位序列相似引起的HiC噪音的必要输入。

构建的方法有很多种，这里列举两个方法。

## 方法一: 基于BLAST

基于BLAST的方法需要你提供一个染色体级别组装的近缘物种。需要提供4个文件，近缘物种的CDS和GFF文件，多倍体物种的CDS和GFF文件。最重要的是，你得保证CDS序列中的编号和GFF文件的基因编号一致，这也是后续处理至关重要的一步。

**问题一**: 我们如何快速获取待组装物种的CDS和GFF文件？最快速的方法就是利用转录组序列基于STAR + StringTie的组合进行拼接，并利用gffread从中提取出cdna序列。

假设你的多倍体文件名为`contig.fa`, 转录组数据为`read_R1.fq.gz`, `read_R2.fq.gz`

```bash
mkdir -p STAR
ref=contig.fa
THREADS=40
fq1=read_R1.fq.gz
fq2=read_R2.fq.gz
prefix=sample
# 建立索引
STAR \
  --runThreadN $THREADS \
  --runMode genomeGenerate \
  --genomeDir STAR \
  --genomeFastaFiles $ref
#比对
STAR \
	--genomeDir STAR \
	--runThreadN $THREADS \
	--readFilesIn $fq1 $fq2 \
	--readFilesCommand zcat \
	--outFileNamePrefix $prefix\
	--outSAMtype BAM SortedByCoordinate \
	--outBAMsortingThreadN $THREADS \
	--outSAMstrandField intronMotif \
	--outFilterIntronMotifs RemoveNoncanonical
# 拼接转录本
stringtie ${prefix}Aligned.sortedByCoord.out.bam -p $THREADS -o qry.gtf
# GTF 转成 GFF3
gffread qry.gtf -o qry.gff
# 提取cdna
gffread qry.gff -g $ref -w qry.fa
```

**问题二**:  后续的处理脚本要求输入序列的编号和GFF的基因编号相同，并且编号在GFF文件中的格式为`NAME=xxx`。

对于我们上一步处理的qry.gff和qry.fa，我们只需要简单粗暴地替换qry.gff里的信息即可。

```bash
sed -e 's/transcript/gene/' -e 's/ID/Name/' qry.gff > qry_gene.gff
```

对于近缘物种的cds序列和gff文件，考虑到后续的脚本只认gene所在的行，而cds里面存放的是转录本序列，我们可以提取mRNA所在行，然后将其改名成gene。假设参考序列是`ref.fa`,  参考的GFF文件名是`ref.gff`

```bash
ref=ref.gff
grep '[[:blank:]]mRNA[[:blank:]]' $ref | sed -e 's/mRNA/gene/' -e 's/ID/Name/' > ref_gene.gff
```

综上你得到了如下四个文件, ref.fa, ref_gene.gff, qry.fa, qry_gene.gff

下面的操作，参考自[官方教程]( https://github.com/tangerzhang/ALLHiC/wiki/ALLHiC:-identify-allelic-contigs )完成，

**第一步**： 建立BLAST索引

```bash
makeblastdb -in ref.fa -dbtype nucl
```

**第二步**: 将我们物种序列比对到近缘物种

```bash
blastn -query qry.fa -db ref.fa -out qry_vs_ref.blast.out -evalue 0.001 -outfmt 6 -num_threads 4 -num_alignments 1
```

**第三步**: 过滤低质量的HIT，此处标准是相似度低于60%, 覆盖度低于80%

```bash
blastn_parse.pl -i qry_vs_ref.blast.out -o Eblast.out -q qry.fa -b 1 -c 0.6 -d 0.8 
```

**第四步**: 基于BLAST结果鉴定等位序列

```bash
classify.pl -i Eblast.out -p 4 -r ref_gene.gff -g qry_gene.gff
# -p 是等位基因的数目
```

最终输出的Allele.ctg.table就是后续ALLHiC的prune输入。

## 方法二:  基于GMAP

第二个方法比较简单，只需要提供多倍体基因组的基因组序列和近缘物种的cds序列即可

**第一步**: 运行GMAP获取GFF3文件

```bash
qry=target.genome  # 你的多倍体序列
ref_cds=reference.cds.fasta # 同源物种的cds
N=4 #你的基因组倍性, 4表示的是4倍体
gmap_build -D . -d DB $qry
gmap -D . -d DB -t 12 -f 2 -n $N  $ref_cds > gmap.gff3
```

**第二步**: 获取  allelic.ctg.table 

```bash
perl gmap2AlleleTable.pl referenece.gff3
```

 gmap2AlleleTable.pl脚本是ALLHiC的一部分。

## 参考资料

- <https://github.com/tangerzhang/ALLHiC/issues/16> 
- <https://github.com/tangerzhang/ALLHiC/wiki/ALLHiC:-identify-allelic-contigs>