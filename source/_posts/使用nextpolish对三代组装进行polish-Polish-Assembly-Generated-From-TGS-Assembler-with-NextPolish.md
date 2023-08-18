---
title: 使用nextpolish对三代组装进行polish
date: 2019-10-25 11:25:33.711
updated: 2019-10-25 11:25:33.711
url: /archives/Polish-Assembly-Generated-From-TGS-Assembler-with-NextPolish
categories: 生信软件工具箱
tags: 组装
---



NextPolish是武汉未来组开发的一个三代基因组polish工具（另外一个常用软件是Pilon）。NextPolish可以使用二代短读序列或者三代序列或者两者结合去纠正三代长读长序列在组装时导致的碱基错误(SNV/Indel)。由于它是专为polish设计，因此在运行速度和内存使用上都优与Pilon。

## 软件安装

先确保自己的服务器上安装了Python2.7， 且有Shutil和Signal，或者你可以利用conda新建一个python2.7的环境。

```bash
# shell
python -V
Python 2.7.15
# Python 交互命令行
import shutil
import signal
```

之后到<https://github.com/Nextomics/NextPolish/releases>找最新的版本，写这篇时的最新版本就是1.05

```bash
mkdir -p ~/opt/biosoft
cd ~/opt/biosoft
wget https://github.com/Nextomics/NextPolish/releases/download/v1.0.5/NextPolish.tgz
tar -zxvf NextPolish.tgz
# 编译软件
cd NextPolish && make -j 10
# 加入到.bashrc或.zshrc
export PATH=~/opt/biosoft/NextPolish:$PATH
```

## 软件使用

**注意**：如果你的基因组用的是miniasm这类缺少consensus步骤的组装软件，那么你需要先用运行如下命令，或者是运行`racon`利用三代序列进行polish。否则，由于基因组上存在过高的错误率，导致二代序列错误比对，影响polish效果。

```bash
threads=20  
genome=input.genome.fa # 组装的基因组
lgsreads=input.lgs.reads.fq.gz # 三代长度序列
# 将三代回帖到参考基因组
minimap2 -a -t ${threads} -x map-ont/map-pb ${genome} ${lgsreads}| \
    samtools view -F 0x4 -b - | \
    samtools sort - -m 2g -@ ${threads} -o genome.lgs.bam  
#建立索引
samtools index -@ ${threads} genome.lgs.bam
samtools faidx ${genome}
# 使用nextPolish.py 进行polish
python ~/opt/biosoft/NextPolish/lib/nextPolish.py \
    -g ${genome} -t 5 --bam_lgs genome.lgs.bam -p ${threads} > genome.lgspolish.fa
```

生成的genome.lgspolish.fa才能用于后续的二代polish步骤。

NextPolish要求我们需要准备两个文件：

- `run.cfg`: 配置文件，
- `sgs.fofn`: 二代测序文件的位置信息

以[使用NextDenovo组装Nanopore数据](/archives/Assembly-nanopore-with-NextDenovo)文章组装的结果为例进行介绍。在分析目录下有三个文件。

- 三代组装结果: nextgraph.assembly.contig.fasta
- 二代序列: ERR2173372_1.fastq,ERR2173372_2.fastq

第一步：创建一个文件，用于记录二代序列的位置信息

```bash
realpath ERR2173372_1.fastq ERR2173372_2.fastq  > sgs.fofn
```

第二步：配置run.cfg文件

```bash
# 从NextPolish目录下复制配置文件
cp ~/opt/biosoft/NextPolish/doc/run.cfg run2.cfg
```

修改配置文件

```bash
[General]
job_type = local
job_prefix = nextPolish
task = default
rewrite = 1212
rerun = 3
parallel_jobs = 2
multithread_jobs = 10
genome = ./nextgraph.assembly.contig.fasta
genome_size = auto
workdir = ./01_rundir
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100
```

其中需要修改的参数为，其余参数查看[官方的参数配置说明](https://github.com/Nextomics/NextPolish/blob/master/doc/OPTION.md):

- `job_type`: 任务类型，local表示单个节点运行。由于NextPolish使用DRMAA进行任务投递，因此还支持，SGE, PBS和SLURM
- `task`: 任务类型， 用12,1212,121212,12121212来设置polish的轮数，建议迭代2轮就可以了。
- `parallel_jobs`和`multithread_jobs`表示同时投递的任务数和每个任务的线程数，此处2 X 10=20
- `genome`: 表示组装基因组的位置
- `workdir`: 输出文件所在目录
- `sgs_options`: 该选项设置二代测序polish的参数，包括-use_duplicate_reads, -unpaired, -max_depth, -bwa, -minimap2(默认使用)

运行方法

```bash
nextPolish run2.cfg &
```

在最后输出日志中，会提示最终存放的文件在什么位置，然后将这些文件合并到单个文件即可。

