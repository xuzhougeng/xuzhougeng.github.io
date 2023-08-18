---
title: 「三代组装」使用新版Falcon进行三代测序基因组组装
date: 2019-08-20 15:53:04.991
updated: 2019-08-24 18:00:46.82
url: /archives/Falcon-genome-assembly-tools
categories: 生信软件工具箱
tags: 组装
---

> 这里的新版指的是PacBio公司在2018年9月发布pb-assembly, 而这篇文章是在2018年9月30日发的。

今年早些时候在参加三代培训时，听说PacBio会在今年对Falcon进行一些改变。前几天我在读 **readthedocs**上的Falcon文档时，发现了文档页面上方出现了这样两栏提醒

![注意](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/attentions-d2f8c343963f471f8a2fcd85a64e708e.png)

其中[pb_assembly](https://github.com/PacificBiosciences/pb-assembly)就是新的FALCON组装套装在GitHub上的使用文档，经过了几天的探索，我对它有了一点了解，写一篇教程作为官方文档的一些补充吧。

## 新版亮点

1. 整合了串联重复序列和遍在重复序列的屏蔽（之前没有这一步）
2. 以GFA格式存放graph文件，后续可以用Bandage进行可视化
3. 通过算法和性能优化，提高了Associate Contigs的准确性
4. 分析流程的性能优化

## 软件安装和数据准备

Falcon终于拥抱了bioconda, 这也就意味着我们再也不需要用到他们原本笨拙的安装脚本，浪费时间在安装软件上。

```bash
conda create -n pb-assembly  pb-assembly
source activate pb-assembly
# 或者
conda create -p ~/opt/biosoft/pb-assembly pb-assembly
source activate ~/opt/biosoft/pb-assembly 
```

这里使用<https://pb-falcon.readthedocs.io/en/latest/tutorial.html>上的所用的E. coli数据集

```bash
wget https://downloads.pacbcloud.com/public/data/git-sym/ecoli.m140913_050931_42139_c100713652400000001823152404301535_s1_p0.subreads.tar.gz
tar -xvzf ecoli.m140913_050931_42139_c100713652400000001823152404301535_s1_p0.subreads.tar.gz
```

解压缩后的文件夹里有三个300M的Fasta文件, 将他们的实际路径记录到`input.fofn`中
```bash
ecoli.1.fasta
ecoli.2.fasta
ecoli.3.fasta
```

## 准备配置文件

为了进行组装，需要准备一个配置文件。我的配置文件为`fc_run.cfg`，内容如下。你们可以先预览一下，后面看我的解释说明。

```bash
#### Input
[General]
input_fofn=input.fofn
input_type=raw
pa_DBdust_option=
pa_fasta_filter_option=pass
target=assembly
skip_checks=False
LA4Falcon_preload=false

#### Data Partitioning
pa_DBsplit_option=-x500 -s50
ovlp_DBsplit_option=-x500 -s50

#### Repeat Masking
pa_HPCTANmask_option=
pa_REPmask_code=1,100;2,80;3,60

####Pre-assembly
genome_size=0
seed_coverage=20
length_cutoff=3000    
pa_HPCdaligner_option=-v -B128 -M24
pa_daligner_option=-e.7 -l1000 -k18 -h80 -w8 -s100
falcon_sense_option=--output-multi --min-idt 0.70 --min-cov 2 --max-n-read 800
falcon_sense_greedy=False

####Pread overlapping
ovlp_daligner_option=-e.96 -l2000 -k24 -h1024 -w6 -s100
ovlp_HPCdaligner_option=-v -B128 -M24

####Final Assembly
overlap_filtering_setting=--max-diff 100 --max-cov 300 --min-cov 2
fc_ovlp_to_graph_option=
length_cutoff_pr=2000

[job.defaults]
job_type=local
pwatcher_type=blocking
JOB_QUEUE=default
MB=32768
NPROC=6
njobs=32
submit = /bin/bash -c "${JOB_SCRIPT}" > "${JOB_STDOUT}" 2> "${JOB_STDERR}"

[job.step.da]
NPROC=4
MB=32768
njobs=32
[job.step.la]
NPROC=4
MB=32768
njobs=32
[job.step.cns]
NPROC=8
MB=65536
njobs=5
[job.step.pla]
NPROC=4
MB=32768
njobs=4
[job.step.asm]
NPROC=24
MB=196608
njobs=1
```

根据注释信息，文件分为"input", "Data partitioning", "Repeat Masking", "Pre-assembly", "Pread overlapping", "Final Assembly", 以及最后的任务调度部分，让我们分别看下这里面的内容

### 输入(Input)

```bash
#### Input
[General]
input_fofn=input.fofn
input_type=raw
pa_DBdust_option=
pa_fasta_filter_option=pass
target=assembly
skip_checks=False
LA4Falcon_preload=false
```

输入这里参数比较简单，基本不需要做任何改动，除了 pa\_fasta_filter_options，用于处理一个ZMW（测序翻译孔）有多条subread时，到底选择哪一条的问题。

- "pass": 不做过滤，全部要。
- "streamed-median": 表示选择大于中位数的subread
- "streamed-internal-median": 当一个ZMW里的subread低于3条时选择最长，多于单条则选择大于中位数的subread

> 0.01版本pb-assembly的`pa_DBdust_option`有一个bug，也就是里面的参数不会传递给DBdust, DBdust是对read进行soft-masking，一般都用默认参数，因此这个bug问题不大。

### 数据分配(Data Partitioning)

```bash
pa_DBsplit_option=-x500 -s50
ovlp_DBsplit_option=-x500 -s50
```

这部分的设置会将参数传递给`DBsplit`，将数据进行拆分多个block，后续的并行计算都基于block。`-s 50`表示每个block大小为50M。 这适用于基因组比较小的物种，如果是大基因组则应该设置为`-s 200`或者`-s 400`

### 重复屏蔽(Repeat Masking)

```bash
#### Repeat Masking
pa_HPCTANmask_option=
pa_REPmask_code=1,100;2,80;3,60
```

屏蔽重复序列可以在不损失组装准确性的同时，提高后续组装的overlap/daligner步骤10~20倍速度，见[Detecting and Masking Repeats](https://dazzlerblog.wordpress.com/2016/04/01/detecting-and-soft-masking-repeats/).

pa\_HPCTANmask_option的参数会传给串联重复步骤的`HPCTANmask`, 而`pa_REPmask_code`很复杂，它分为三次迭代，因此这里1:100;2,80;3,60 就表示第一次迭代检测每个block中出现超过100次的序列，第二次迭代将2个block合并一起检测超过80次的序列，第三次将3个block进行合并检测超过60次的序列。

> <https://github.com/PacificBiosciences/pbbioconda/issues/20>

### 预组装（纠错）pre-assembly

```bash
####Pre-assembly
genome_size=0
seed_coverage=20
length_cutoff=3000    
pa_HPCdaligner_option=-v -B128 -M24
pa_daligner_option=-e.7 -l1000 -k18 -h80 -w8 -s100
falcon_sense_option=--output-multi --min-idt 0.70 --min-cov 2 --max-n-read 800
falcon_sense_greedy=False
```

当`length_cutoff=-1`时，设置`genome_size`和`seed_coverage`会自动计算要过滤的序列。否则是过滤低于一定长度的read。

pa\_HPCdalinger_option参数不需要调整，-M表示每个进程的内存为24G，一般200M的block对应16G。

pa\_daligner_option的参数比较重要：

- `-e`:错误率，低质量序列设置为0.70，高质量设置为0.80。 值越高避免单倍型的坍缩
- `-l`:  最低overlap的长度，文库比较短时为1000， 文库比较长为5000.
- `-k`:  低质量数据为14，高质量数据为18
- `-h`: 表示完全match的k-mer所覆盖的碱基数。和`-l`, `-e`有关，越大越严格。预组装时最大也不要超过最低overlap长度的1/4. 最低就设置为80

### 纠错后相互比对

```bash
####Pread overlapping
ovlp_daligner_option=-e.96 -l2000 -k24 -h1024 -w6 -s100
ovlp_HPCdaligner_option=-v -B128 -M24
```

和上面的参数类似，但是-e的范围调整为0.93-0.96，-l范围调整为1800-6000， -k调整为 18-24

### 最后组装

```bash
overlap_filtering_setting=--max-diff 100 --max-cov 300 --min-cov 2
fc_ovlp_to_graph_option=
length_cutoff_pr=2000
```

这里的参数就可以随便调整了，因为这一步速度很快。 例如`length_cutoff_pr`就可以从2000，提高到15000.

最后还有一部分是任务投递系统，如果是单节点运行，需要注意设置 njobs，这是同时投递的任务数。假如你将[job.step.cns]按照如下的方式设置，那么同时会出现 $ 8 \times 50 = 400 $ 个任务，如果你的内存只有128G，运行一段时间后你的所有内存就会被耗尽，那么基本上你就只能重启服务器了。

```bash
[job.step.cns]
NPROC=8
MB=65536
njobs=50
```

## 运行结果

用上述的配置文件，以`fc_run fc_run.cfg`运行后，最后的`2-asm-falcon/p_ctg.fa`的序列数有4条，最长为4,685,024, 之后我将`length_cutoff_pr`调整为15k，`2-asm-falcon/p_ctg.fa`序列只有一条，长度为`4,638,411`

下载`asm.gfa`到本地，用[Bandage](https://rrwick.github.io/Bandage/)可视化，可以发现组装效果不错。

![Bandage](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/bandage-2c760dfa117940389cedc639a8aa9906.png)

