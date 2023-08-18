---
title: 使用Pilon对基因组进行polish
date: 2019-08-27 20:21:45.5
updated: 2019-08-27 20:23:15.374
url: /archives/Polish-genome-assembly-by-pilon
categories: 生信软件工具箱
tags: 
---

# 使用Pilon对基因组进行polish

## 软件安装

官方提供了编译好的jar包，方便使用

```bash
wget https://github.com/broadinstitute/pilon/releases/download/v1.22/pilon-1.22.jar
java -Xmx16G -jar pilon-1.22.jar
```

如果要顺利运行程序，要求JAVA > 1.7, 以及根据基因组大小而定的内存，一般而言是1M大小的基因对应1GB的内存。

## 总览

Pilon有如下作用

1. 对初步组装进行polish
1. 寻找同一物种不同株系间的变异，包括结构变异检测

他以FASTA和BAM文件作为输入，根据比对结果对输入的参考基因组进行提高，包括

- 单碱基差异
- 小的插入缺失(indels)
- 较大的插入缺失或者block替换时间
- 填充参考序列中的N
- 找到局部的错误组装

最后它输出polish后的FASTA文件, 以及包含变异信息的VCF文件(可选)

## 分析流程

推荐使用PCR-free建库测序得到的Illumina paired-end数据，这样子避免了PCR-duplication,有效数据更多，也不需要在分析过程中标记重复。

下面步骤，假设你的组装文件为`draft.fa`, 质量控制后的illumina双端测序数据分别为`read_1.fq.gz`和`read_2.fq.gz`

第一步：比对

```bash
bwa index -p index/draft draft.fa
bwa mem -t 20 index/draft read_1.fq.gz read_2.fq.gz | samtools sort -@ 10 -O bam -o align.bam
samtools index -@ 10 align.bam
```

第二步：标记重复（非PCR-free建库)

```bash
sambamba markdup -t 10 align.bam align_markdup.bam
```

第三步：过滤高质量比对的read

```bash
samtools view -@ 10 -q 30 align_markdup.bam > align_filter.bam
samtools index -@ 10 align_filter.bam
```

第三步：使用Pilon进行polish

```bash
MEMORY= #根据基因组大小而定
java -Xmx${MEMORY}G -jar pilon-1.22.jar --genome draft.fa --frags align_filer.bam \
    --fix snps,indels \
    --output pilon_polished --vcf &> pilon.log
```

关于Pilon的一些参数说明：

- `--frags`表示输入的是1kb以内的paired-end文库，`--jumps`表示 大于1k以上的mate pair文库,  `--bam`则是让软件自己猜测
- `-vcf`: 输出一个vcf文件，包含每个碱基的信息
- `--fix`:  Pilon将会处理的内容，基本上选`snps`和`indels`就够了
- `--variant`: 启发式的变异检测，等价于`--vcf --fix all,breaks`, 如果是polish不要使用该选项
- `minmq`: 用于Pilon堆叠的read最低比对质量，默认是0。

## 阅读日志输出

> 这个日志文件是标准输出而不是标准错误输出，不过保险起见用`&>`

最开始，Pilon会输出他的版本信息（如下示例），以及将会对基因组做的调整,

```bash
Pilon version 1.14 Sat Oct 31 14:30:00 2015 -0400
Genome: genome.fasta
Fixing snps, indels
```

其中Fixing后面的含义为：

- "snps":  单碱基差异
- "indels":小的indel的差异
- "amb":  替换原有的N
- "gaps": 填充基因组的gap
- "local"： 检测和修改错误组装
- "all":  上述所有
- "none": 不是上述的任何一种

接着Pilon会分染色体对BAM文件进行处理，根据BAM文件进行堆叠(pileup), 这个时候它会输出有效reads的深度，这里的有效reads包括未成对的read和正确成对的read。

```bash
Processing ctg1:1-5414473
frags align_mkdup.bam: coverage 19
Total Reads: 808985, Coverage: 19, minDepth: 5
```

从Pilon v1.4开始，Pilon还会输出基因组得到确认的碱基比例。

```bash
Confirmed 5403864 of 5414473 bases (99.80%)
```

后续是Pilon将会对原参考基因组做的一些调整的总体情况，如下表示纠正2个snp, 2个小的插入，4个缺失。

```bash
Corrected 2 snps; 0 ambiguous bases; corrected 2 small insertions totaling 12 bases, 4 small deletions totaling 6 bases
```

最后声明当前部分处理结束

```bash
Finished processing ctg1:1-5414473
```

如果，在`--fix`中选了`gaps`, 那么输出的内容还有如下内容。其中`82048 -0 +276`解释为在坐标82428处移除0个碱基，插入276个碱基。

```bash
# Attempting to fill gaps
fix gap: scaffold00001:82428-93547 82428 -0 +276 ClosedGap


```

参考资料

- <https://github.com/broadinstitute/pilon/wiki/Standard-Output>
- <https://github.com/broadinstitute/pilon/wiki/Methods-of-Operation>