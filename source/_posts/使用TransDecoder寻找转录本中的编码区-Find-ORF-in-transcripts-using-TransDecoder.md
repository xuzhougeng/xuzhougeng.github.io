---
title: 使用TransDecoder寻找转录本中的编码区
date: 2019-08-28 22:43:45.469
updated: 2019-09-16 15:05:12.362
url: /archives/Find-ORF-in-transcripts-using-TransDecoder
categories: 生信软件工具箱
tags: 组装 | 转录组 | 注释
---


TransDecoder能够从转录本序列中鉴定候选编码区。这些转录本序列可以来自于Trinity的从头组装，或者来自于Cufflinks或者StringTie的组装结果。

## 软件安装

从<https://github.com/TransDecoder/TransDecoder/releases>下载最新版的TransDecoder，以v5.5.0为例

```bash
mkdir -p ~/opt/biosoft && cd ~/opt/biosoft
wget https://github.com/TransDecoder/TransDecoder/archive/TransDecoder-v5.5.0.zip
unzip TransDecoder-v5.5.0.zip
mv TransDecoder-TransDecoder-v5.5.0 TransDecoder-v5.5.0
```

## 运行TransDecoder

我们从cufflinks或stringtie输出的gtf文件开始分析流程，因此你会有两个输入文件

- transcripts.gtf: 记录预测转录本的GTF文件
- genome.fasta: 参考基因组序列

**第一步**: 从GTF文件中提取FASTA序列

```bash
~/opt/biosoft/TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl transcripts.gtf genome.fasta > transcripts.fasta
```

**第二步**: 将GTF文件转成GFF3格式

```bash
~/opt/biosoft/TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl transcripts.gtf > transcripts.gff3
```

**第三步**: 预测转录本中长的开放阅读框, 默认是100个氨基酸，可以用`-m`修改

```bahs
~/opt/biosoft/TransDecoder-v5.5.0/TransDecoder.LongOrfs -t transcripts.fasta
```

**第四步**: 使用DIAMOND对上一步输出的`transcripts.fasta.transdecoder.pep`在蛋白数据库中进行搜索，寻找同源证据支持

```bash
# 下载数据并解压缩
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
# 建立索引
diamond makedb --in uniprot_sprot.fasta --db uniprot_sprot.fasta
# BLASTP比对
diamond blastp -d uniprot_sprot.fasta -q transcripts.fasta.transdecoder_dir/longest_orfs.pep --evalue 1e-5 --max-target-seqs 1 > blastp.outfmt6
```

> 关于DIAMOND的使用，参考这篇说明[DIAMOND: 超快的蛋白序列比对软件](/archives/Fast-and-sensitive-protein-alignment-using-diamond)

**第五步**: 预测可能的编码区

```bash
~/opt/biosoft/TransDecoder-v5.5.0/TransDecoder.Predict \
    -t transcripts.fasta \
    --retain_blastp_hits blastp.outfmt6 
```

**第六步**: 生成基于参考基因组的编码区注释文件

```bash
~/opt/biosoft/TransDecoder-v5.5.0/util/cdna_alignment_orf_to_genome_orf.pl \
     transcripts.fasta.transdecoder.gff3 \
     transcripts.gff3 \
     transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3
```

最终输出文件如下:

- transcripts.fasta.transdecoder.pep: 最终预测的CDS对应的蛋白序列
- transcripts.fasta.transdecoder.cds: 最终预测的CDS序列
- transcripts.fasta.transdecoder.gff3: 最终ORF对应的GFF3
- transcripts.fasta.transdecoder.bed: 以BED格式存放ORF位置信息
- transcripts.fasta.transdecoder.genome.gff3: 基于参考基因组的GFF3文件

其中BED和GFF3可以放到IGV上展示，手动检查下结果

> 假如是Trinity从头预测的转录本，没有参考基因组，那么就运行第三步，第四步和第五步

在transcripts.fasta.transdecoder.cds文件中，每个fasta序列的头信息部分会有一个ORF type，分为如下几个类型

- 'complete': 包含起始密码子和终止密码子
- '5prime_partial': 缺失起始密码子, 可能只有部分的N端序列
- '3prime_partial': 缺失终止密码子, 可能只有部分的C端序列
- 'internal': 意味着同时是5prime-partial和3prime-partial

通常而言，我们会用complete的cds寻找同源证据，然后选择高可信度的序列用于训练，而不是哪些特别长且没有已知蛋白支持的序列。

## 参考资料

- https://github.com/TransDecoder/TransDecoder/wiki
- https://github.com/PASApipeline/PASApipeline/wiki/PASA_abinitio_training_sets