---
title: 使用BRAKER2进行基因组注释
date: 2019-10-10 18:57:19.681
updated: 2020-06-28 07:21:49.532
url: /archives/genome-annotation-with-braker2
categories: 生信软件工具箱
tags: 注释 | 流程工具
---

# 使用BRAKER2进行基因组注释

BRAKER2是一个基因组注释流程，能够组合GeneMark，AUGUSTUS和转录组数据。

在使用软件之前，有几点需要注意下

- 尽量提供高质量的基因组。目前随着三代测序价格下降，这一点问题不大。
- 基因组命名应该简单，最好就是">contig1"或">tig000001"
- 基因组需要屏蔽重复序列
- 默认参数通常表现效果就很好，但是也要根据物种来
- 一定要对注释结果进行检查，别直接使用

## 软件安装

BRAKER的依赖软件不少，且Perl需要安装的模块也很多，但是我们可以直接用bioconda进行安装

```bash
conda create -n braker2 braker2
```

使用conda安装结束后会有一些提示语，总结如下

- 保证AUGUSTUS的config目录能够有可写权限（自己用conda安装不需要考虑这个问题）
- GeneMark/GenomeThreader/ProtHint需要额外下载安装

> 尽管可以使用容器化技术，但是由于软件运行时需要对AUGUSTUS里的配置文件进行读写，需要额外设置，因此不在此介绍

由于部分软件的版权限制，也就是GeneMakr, GenomeThreader, ProHint，我们还需要手动安装这些软件。但其中只有GeneMark是必须的，因为我们更多是利用RNA-seq数据进行模型训练，而GenomeThreader, ProHint是利用同源蛋白进行注释才用到，其中GenomeTrheads是处理近源蛋白，而ProHint是处理未知距离的蛋白。

从 <http://exon.gatech.edu/GeneMark/license_download.cgi> 下载GeneMark，并按照文档进行安装，最后添加环境变量

```bash
export GENEMARK_PATH=/your_path_to_GeneMark-ET/gmes_petap/
#例如,我的安装路径为/opt/biosoft/gm_et_linux_64
export GENEMARK_PATH=/opt/biosoft/gm_et_linux_64
```

安装完成之后，建议运行下面这一步检查软件依赖

```bash
braker.pl --checkSoftware
```

## 软件运行

BRAKER根据数据类型，有[不同的运行模式](https://github.com/Gaius-Augustus/BRAKER#overview-of-modes-for-running-braker)，但根据现状其实最常见的情况是测了一个基因组，并且还测了二代的转录组。因此假设你手头有下面这些数据

- 基因组序列: genome.fasta
- 转录组数据: XX_R1.fq.gz, XX_R2.fq.gz, YY_R1.fq.gz, YY_R2.fq.gz, ZZ_R1.fq.gz, ZZ_R2.fq.gz

第一步: 屏蔽基因组中的重复序列，这一步参考[使用RepeatModeler和RepeatMasker注释基因组重复序列](/archives/Annotate-repeat-in-genome-with-RepeatModeler-and-RepeatMasker)

```bash
RepeatMasker -xsmall -species arabidopsis -pa 40 -e ncbi  -dir . genome.fasta
#-xsmall: soft-mask
```

这一步输出的genome.fasta.masked将是后续注释的输入

第二步: 使用STAR将FastQ比对到参考基因组，STAR使用说明参考[「RNA-seq分析软件」RNA-seq比对工具STAR学习笔记](/archives/RNA-seq-aligner-STAR)

```bash
mkdir -p STAR
# 建立索引
STAR \
    --runThreadN 20 \
    --runMode genomeGenerate \
    --genomeDir STAR \
    --genomeFastaFiles genome.fasta
# 比对
STAR \
	--genomeDir STAR \
	--runThreadN 20 \
	--readFilesIn XX_R1.fq.gz, XX_R2.fq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix xx_ \
	--outSAMtype BAM SortedByCoordinate \
	--outBAMsortingThreadN 10 \
	--outSAMstrandField intronMotif \
	--outFilterIntronMotifs RemoveNoncanonical
mv xx_Aligned.sortedByCoord.out.bam xx.bam
```

输入结果为 xx.bam 如果测了多个组装的转录组，为每个样本运行一次比对生成多个BAM文件。

第三步: 运行BRAKER2。`--cores 40`表示使用40个线程, `--softmaksing`表示输入序列是软屏蔽。

```bash
braker.pl --cores 40 --species=yourSpecies --genome=genome.fasta.masked \
     --softmasking --bam=xx.bam,yy.bam,zz.bam \
     --gff3
```

这里xx.bam, yy.bam,zz.bam指的是不同组织的BAM文件，需要根据实际情况进行修改。

对于链特异性转录组测序结果，需要将BAM文件拆分成正链和负链两个bam，然后设置`--stranded`参数。

```bash
braker.pl --species=yourSpecies --genome=genome.fasta.masked \
   --softmasking --bam=plus.bam,minus.bam,unstranded.bam \
   --stranded=+,-,. --UTR=on
```

如果需要使用近源序列作为输入(protein.fa)，需要加上`--prot_seq`和`--prg`参数，速度会慢一些

```bash
braker.pl --cores 40 --species=yourSpecies --genome=genome.fasta.masked \
     --softmasking --bam=xx.bam,yy.bam,zz.bam \
     --gff3 \
     --prot_seq=proteins.fa --prg=exonerate \
```

> 个人感觉，只要转录组测得够，不需要考虑同源蛋白数据。

最终会输出蛋白序列和CDS序列以及GFF文件

- augustus.hints.gtf: 在`--esmode`时不会出现
- augustus.hints_utr.gtf: 在设置`--UTR=on `时生成
- augustus.ab_initio.gtf: 在设置`--AUGUSTUS_ab_initio`生成
- augustus.ab_initio_utr.gtf: 在从头预测的基础上设置`--UTR=pn`时生成
- GeneMark-E*/genemark.gtf: GeneMark的注释结果
- braker.gtf: AUGUSTUS或GeneMark预测的所有结果
- hintsfile.gff: 来自于RNA-seq/蛋白数据的外部证据信息

此外，如果没有设置`--skipGetAnnoFromFasta`, 还会有`augustus.hints.aa`和`augustus.hints.codingseq`

**注意**: BRAKER没有断点重跑功能，因此如果中途运行中断，重新运行之前的命令并不会断点重跑，反而会因为同样的参数和之前运行在AUGUSTUS的config/species生成的目录冲突而中断。另外设置`--useexisting`只是覆盖之前`config/species`对应的目录，而不是利用已有的文件重新开始

尽管如此，我们还是能够在BRAKER中通过跳过一些步骤来避免一些重复运算

- `--skipGeneMark-ES/--skipGeneMark-ET/--skipGeneMark-EP/--skipGeneMark-ETP`: 跳过GeneMark训练这一步，使用之前的结果`braker/GeneMark-ET/genemark.gtf`，或者使用`--geneMarkGtf`设置GTF文件
- `--skipAllTraining`: 跳过GeneMark和AUGUSTUS的模型训练，使用之前的参数和文件(`--useexisting`)运行AUGUSTUS

官方教程也提供了一些案例，见<https://github.com/Gaius-Augustus/BRAKER#starting-braker-on-the-basis-of-previously-existing-braker-runs>

除此之外，如果想要提高速度，还可以设置如下参数

- `--skipOptimize`:  跳过参数优化步骤（不推荐）
- `--rounds`: 参数优化迭代数，默认是5，追求速度可以设置的低一些

## 可能问题

运行时出现和OpenBLAS相关的警告和报错

```bash
OpenBLAS blas_thread_init: RLIMIT_NPROC 4096
OpenBLAS blas_thread_init: pthread_create failed for thread 26 of 128: Resource temporarily unavailable
```

需要在运行前设置如下两个环境变量

```bash
export OMP_NUM_THREADS=1 #限制CPU数
export USE_SIMPLE_THREADED_LEVEL3=1
```

使用conda安装时可能会出现的问题

```bash
Error in file bamToWig.py at line 172: Return code of subprocess was 127
```

原因是因为`faToTwoBit`程序出错

```bash
faToTwoBit
faToTwoBit: error while loading shared libraries: libssl.so.1.0.0: cannot open shared object file: No such file or directory
```

这是因为conda没能正确处理依赖关系，openssl版本过高，解决方法如下

```bash
# 建立软链接
cd ~/miniconda3/envs/braker2/lib
ln -s libssl.so libssl.so.1.0.0
ln -s libcrypto.so libcrypto.so.1.0.0
```

## 参考资料

- BRAKER2官方教程: https://github.com/Gaius-Augustus/BRAKER