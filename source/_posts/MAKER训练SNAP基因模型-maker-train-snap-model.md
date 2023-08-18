---
title: MAKER训练SNAP基因模型
date: 2020-07-06 00:09:21.893
updated: 2020-07-07 07:15:54.553
url: /archives/maker-train-snap-model
categories: 基因组学
tags: 注释 | MAKER
---


## 准备阶段

训练SNAP模型，需要准备三个文件，分别是参考基因组序列，组装的转录本序列和同源蛋白序列。

对于参考基因组序列，我们要保证N50需要超过预期基因长度的中位数，否则注释效果不好。不过目前的基因组在三代测序的加持下，基本上都不是问题。

对于同源蛋白， 建议只用UniProt/Sprot的人工检查过的高质量蛋白序列，而不是盲目参考文献，使用近源物种的所有蛋白。除非你的近源物种都是模式物种，蛋白序列可靠性高，否则用错误的输入进行训练，数据越多反而错的越多。

我们可以在<ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/>选择合适的uniprot_sprot数据,  然后将其输出为fasta格式。以植物为例

```bash
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_plants.dat.gz
zcat uniprot_sprot_plants.dat.gz |\
    awk '{if (/^ /) {gsub(/ /, ""); print} else if (/^AC/) print ">" $2}' |\
    sed 's/;$//'> protein.fa
```

对于转录本，我们通常会测一些转录组数据，有三种策略可以得到转录本。（这里暂时不考虑三代全长转录本）

1. Trinity重头组装转录本
1. 使用STAR + Trinity 获取转录本
1. 使用STAR + StringTie  + gffread 获取转录本

对于这三种策略，不推荐策略一，因为在有参考基因组的情况下，策略二不但计算效率高，而且能避免组装错误（多倍体等位基因之间相似度高）。对于策略二和策略三，我会推荐策略三。因为对于靠的比较近的基因，Trinity很可能会把这两个基因装成一个。

![Fig1](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/07/image-a2927f2e13ff4a24b7fab8da45207278.png)

并且利用该转录本作为输入训练SNAP模型，之后以SNAP模型作为输入，将转录组和同源蛋白作为证据而不是直接用作模型，我们再检查maker的结果, 也会发现使用StringTie进行组装的结果才是对的。

![Fig2](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/07/image-5197490a14354bdd9e58e75a59468114.png)

因此使用策略二不但计算量大，而且有些情况下还会导致过近的转录本错误融合，反而影响了最终效果，因此我最终推荐**策略三**。当然，这是我定性通过IGV浏览结果得出的结论，样本小，结论未必可靠，仅供参考。

## 训练阶段

假设我们准备的三个文件分别命名为, genome.fa, Trinity-GG.fasta 和 protein.fa

接着使用`maker -CTL`新建配置文件, 设置如下选项

```bash
genome=genome.fa
est=Trinity-GG.fasta
protein=protein.fa
est2genome=1
protein2genome=1
```

然后使用`mpiexec -n 线程数 maker &> run.log`运行程序。

处理结果后，我们新建一个snap目录训练模型

```bash
mkdir snap && cd snap
gff3_merge -d ../genome.maker.output/genome_master_datastore_index.log
```

使用makerzff构建输入文件

```bash
maker2zff -c 0.8 -e 0.8 -o 0.8 -x 0.2 genome.all.gff
```

`maker2zff`会根据AED(-x)和QI值进行过滤，其中QI值一共有9项，每一项的含义如下

1. Length of the 5' UTR
1. Fraction of splice sites confirmed by an EST alignment (`-c`)
1. Fraction of exons that overlap an EST alignmetn(`-e`)
1. Fraction of exons that overlap EST or Protein alignments(`-o`)
1. Fraction of splice site confrimed by a ab-initio prediction(`-a`)
1. Fraction of exons that overlap a ab-initio prediction(`-t`)
1. Number of exons in the mRNA
1. length of the 3' UTR
1. Length of the protein sequence produced by the mRNA (`-l`)

如果QI值第二项为-1，表示没有支持该剪切位点的EST证据.

接着构建模型

```bash
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl snap . > ../snap.hmm
```

修改配置，然后重新运行。MAKER会自动处理冲突的部分，避免重复序列屏蔽等的一些重复计算。

```bash
genome=genome.fa
est=Trinity-GG.fasta
protein=protein.fa
snap=snap.hmm
est2genome=0
protein2genome=0
```

根据输出结果再一次训练模型

```bash
mkdir snap2 && cd snap2
gff3_merge -d ../genome.maker.output/genome_master_datastore_index.log
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl snap . > ../snap2.hmm
```

通常迭代2-3次就够了，毕竟我们可能还会训练AUGUSTUS和GeneMark模型，通过比较多个模型来得到最终结果。

参考资料: <http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018#Training_ab_initio_Gene_Predictors>