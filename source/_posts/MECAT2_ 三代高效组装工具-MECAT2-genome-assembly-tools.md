---
title: MECAT2:三代高效组装工具
date: 2019-08-21 14:16:16.702
updated: 2019-08-24 18:00:16.321
url: /archives/MECAT2-genome-assembly-tools
categories: 生信软件工具箱
tags: 组装
---

# MECAT2: 三代高效组装工具

MECAT2是三代测序数据PacBio的高效组装工具，是之前[MECAT](https://github.com/xiaochuanle/MECAT)的改进版, 修复了之前的很多bug, 使用基于string graph的`fsa`替换了之前的`mecat2canu`。

MECAT2由4个模块组成:

- `mecat2pw`: SMART reads快速和准确地配对比对工具
- `mecat2ref`: SMART reads的参考基因组比对工具
- `mecat2cns`: 基于配对的重叠对存在噪音的reads进行纠错
- `fsa`: 基于string graph的组装工具

目前MECAT2只支持PacBio数据，对于Nanopore数据，肖老师他们开发了`NECAT`

## 安装

MECAT2软件安装非常简单，不存在上一代MECAT的安装难问题了

```r
mkdir -p ~/opt/biosoft
cd ~/opt/biosoft
git clone https://github.com/xiaochuanle/MECAT2.git
cd MECAT2
make -j 4
```

唯一要注意的是尽量要用新版的perl，经测试，Perl 5.26 是可以正常运行脚本。 

最后加入到环境变量`PATH`中

```bash
export PATH=~/opt/biosoft/MECAT2/Linux-amd64/bin:$PATH
```

## 软件使用

我们以 _E. coli_ 数据集作为案例，介绍MECAT2应该如何使用。

**第一步**，我们要下载数据

```bash
mkdir -p ~/ecoli
cd ~/ecoli
wget http://gembox.cbcb.umd.edu/mhap/raw/ecoli_filtered.fastq.gz
gunzip ecoli_filtered.fastq.gz
```

**注意**: 目前MECAT2还不支持gz压缩文件，如果你直接用gz作为输入，会提示core dumped。

MECAT2不是根据文件名来确定输入的序列格式, 也就是如果你把 ecoli_filtered.fastq 命名为  ecoli_filtered.fasta, 组装也不会出错。

**第二步**, 准备模板的config文件

```bash
mecat.pl config ecoli_config_file.txt
```

使用`vim`修改config文件

```bash
PROJECT=ecoli
RAWREADS=ecoli_filtered.fastq
GENOME_SIZE=4800000
THREADS=4
MIN_READ_LENGTH=2000
CNS_OVLP_OPTIONS=""
CNS_OPTIONS="-r 0.6 -a 1000 -c 4 -l 2000"
TRIM_OVLP_OPTIONS="-B"
ASM_OVLP_OPTIONS="-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400"
FSA_OL_FILTER_OPTIONS="--max_overhang=-1 --min_identity=-1"
FSA_ASSEMBLE_OPTIONS=""
USE_GRID=false
CLEANUP=0
CNS_OUTPUT_COVERAGE=30
GRID_NODE=0
```

**第三步**: 原始数据纠错

```bash
mecat.pl correct ecoli_config_file.txt
```

**第四步**: 对纠错后的reads进行组装

```bash
mecat.pl assemble ecoli_config_file.txt
```

输出结果:

- 纠错后的reads： `ecoli/1-consensus/cns_reads.fasta`.
- 最长30X纠错后用于trimming的reads:  `ecoli/1-consensus/cns_final.fasta`.
- trimmed reads:  `ecoli/2-trim_bases/trimReads.fasta`
- 组装的contigs: `ecoli/4-fsa/contigs.fasta`

最终的组装结果为

```bash
Count: 1
Tatal: 4630437
Max: 4630437
Min: 4630437
N25: 4630437
L25: 1
```

能够完美的就装出一条序列

## 配置文件介绍

MECAT2的组装过程的参数都在配置文件中, 也就是之前的`ecoli_config_file.txt`。我按照不同部分对这些参数进行介绍

### 基本参数

这些参数属于改起来不怎么纠结的参数

- `PROJECT=ecoli`: 项目名，之后的所有输出文件都在该项目名的目录下
- `RAWREADS=`: 原始数据的位置，可以是Fasta或者是Fastq， H5格式要先转成Fasta (参考 [Input Format](https://github.com/xiaochuanle/MECAT2#S-input-format)).
- `GENOME_SIZE=`: 基因组大小，单位是bp，如果是120M，那么就是120000000.
- `THREADS=`: 线程数
- `MIN_READ_LENGTH=`: 用于纠错和trim的reads的最低长度. 数据质量好，就长一点，比如说1000到2000
- `USE_GRID=false`: 是否有多个计算节点
- `CLEANUP=0`: 运行结束后删除MECAT2的中间文件, 大基因组的临时文件很大，所以要设置为1.
- `GRID_NODE=0`, 当`USE_GRID=1`时，设置用到的计算节点数，如果是单节点服务器，不需要设置。

### 纠错和修剪阶段

纠错阶段(correct stage)和修剪阶段(trimming stage)，MECAT2调用的`mecat2pw`和`mecat2cns`， 与之相关的配置如下

- `CNS_OVLP_OPTIONS=""`, 在纠错阶段是检测候选重叠的参数, 会传给`mecat2pw`
- `CNS_OPTIONS=""`：原始reads纠错参数，会传递给`mecat2cns`, 
- `TRIM_OVLP_OPTIONS=""`: 在trim阶段，用于检测重叠的参数,会传给`v2asmpm`

但是我发现`v2asmpm`没有文档说明，和软件开发者讨论之后，给出的说明如下

```bash
V2asmpm -Pworkpath -Tt -Sx -Ey -B

-Pworkpath 工作目录是workpath
-Tt 用t个CPU线程
-Sx -Ey 这两个参数一起表示只计算第x到第y卷(包括第y卷), 通常工作目录下会有000001.fasta, 000002.fasta这样的分卷
-B 如果有这个选项, 就输出二进制格式的比对结果; 如果没有这个选项, 就输出文本格式的比对结果
```
因此用默认的`-B`就行了。

### 组装阶段

组装阶段相关的参数如下

- `CNS_OUTPUT_COVERAGE=30`: 选择多少覆盖度的最长纠错后reads进行trim和组装。举个例子，4.8M的 _E. coli_ 30X, 是144MB。 一般选择20~40之间就行， 会传递给`mecat2elr`
- `ASM_OVLP_OPTIONS=""`: 在组装阶段，用于检测重叠的参数，传给`v2asmpm.sh`
- `FSA_OL_FILTER_OPTIONS=""`, 过滤重叠的参数，传递给`fsa_ol_filter`
- `FSA_ASSEMBLE_OPTIONS=""`, 组装trimm reads的参数, 传给`v2asmpm`

这些参数调整起来都很复杂，只能建议看函数帮助文档了。举个例子`FSA_OL_FILTER_OPTIONS="--max_overhang=-1 --min_identity=-1"`根据阅读`fsa_ol_filter`的帮助发现，`-1`表示让程序自己决定

而`ASM_OVLP_OPTIONS="-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400"`， 根据和软件开发者的讨论，是**不需要额外调整**，用默认的即可。

MECAT2是基于比对纠错的组装工具，组装速度上比同类型的Canu和Falcon都快。仅看组装结果的N50参数，MECAT2结果和Canu, Falcon结果都是差不多的。

考虑到BUSCO完整度，Canu是三者最高。不过用三代数据和二代数据进行几轮Polish后会提高到Canu相同水平。三个软件上的BUSCO（embryophyta_odb10）值如下:

```bash
Canu: C:98.8%[S:96.3%,D:2.5%],F:0.5%,M:0.7%,n:1375
MECAT2: C:93.5%[S:91.6%,D:1.9%],F:4.6%,M:1.9%,n:1375
Falcon: C:96.8%[S:95.1%,D:1.7%],F:2.0%,M:1.2%,n:1375
```

因此，如果你原本你是用Falcon做组装，那么可以改成MECAT2了，效果差不多，速度还更快了。对于越大的物种，在有限的资源下，更推荐用MECAT2。

## 参考资料

- <https://github.com/xiaochuanle/MECAT2>