---
title: Oatk：利用HiFi Read进行细胞器组装
date: 2023-06-24 06:24:54.906
updated: 2023-06-24 06:24:54.906
url: /archives/oatk-assembly-of-organelles-using-hifi-read
categories: 基因组学
tags: 基因组
---

Oatk是我见过的一个非常实用、高效的**植物细胞器**基于PacBio HiFi测序的基因组组装工具。这里有两个限定语，一个是植物，一个是PacBio HiFi。如果不满足这两个需求，那就无法使用这个工具。

Oatk包括三个组件，`syncasm`用于从HiFi read中组装细胞器基因组。`hmm_annotation`会利用OatkDB的HMM profile对组装结果进行组装。最后的`pathfinder` 会根据注释信息和组装结果，得到环状的DNA序列。

Oatk软件编译本身仅有一个依赖，zlib，下载源码，就可以编译安装。

```Bash
git clone https://github.com/c-zhou/oatk.git
cd oatk
make -j
```

在实际运行过程中，由于需要用到nhmmscan ，因此还需要额外安装了HMMER，可以使用conda安装。

另外还需要一个OatkDB，里面是作者预处理过的HMM profile，

```Bash
git clone https://github.com/c-zhou/OatkDB.git
```

你可以设置一个环境变量，hmm_db_dir，用于后续分析。

```Bash
hmm_db_dir=/path/to/OatkDB/v20230210/
```

Oatk的运行只需要一个命令，如下

```Bash
input_hifi=/path/to/your/hifi_reads.fq

oatk -k 1001 -c 150 -t 8  \
    -m $hmm_db_dir/angiosperm_mito.fam \
    -p $hmm_db_dir/angiosperm_pltd.fam \
    -o oatk $input_hifi
```

其中-k和-c用于设置syncmer大小和深度。所谓的syncmer，类似Minimizers ，用于从序列中高效的选取K-mers（也就是长度为K的序列）。因此-k 1001 表示，kmer长度最小为1001 bp， -c 150 表示该类型的kmer至少要出现150次。

在后续自己分析设置参数时，-k可以保持1001，-c必须根据实际的测序深度来，作者推荐该值在平均测序深度的5-10倍，用于有效的过滤核基因组序列，确保只筛选出细胞器的序列。

-t 时线程数，根据实际的情况调整。 

-m和-p对应的你实际物种的线粒体和叶绿体的HMM profile。

-o则是输出文件名的前缀。



最终的输出结果为

- oatk.mito.fasta： 线粒体的fasta
- oatk.pltd.fasta：叶绿体的fasta



整体来说，oatk和作者之前开发的工具一样，非常的简洁，高效，并且实用。正如他所用的编程语言C一样。如果你用Hi-C用来搭基因组的话，强烈推荐用他的yahs。



