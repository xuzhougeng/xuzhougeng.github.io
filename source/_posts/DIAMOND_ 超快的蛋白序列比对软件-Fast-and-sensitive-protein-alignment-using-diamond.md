---
title: DIAMOND:超快的蛋白序列比对软件
date: 2019-08-28 21:27:53.717
updated: 2019-09-03 12:21:30.305
url: /archives/Fast-and-sensitive-protein-alignment-using-diamond
categories: 生信软件工具箱
tags: 注释 | 序列比对
---



今天用BLASTX将我的转录本序列在UniProt蛋白数据库(700w条序列)中搜索，80个线程，过了1小时大概就分析1000条吧。实在是有点慢，于是我想到之前耳闻的DIAMOND，据说速度非常快，于是我测试了下。没想到，这工具居然那么快。

根据DIAMOND介绍，它有以下特点

- 比BLAST快500到20,000倍
- 长序列的移框联配分析(frameshift alignment)
- 资源消耗小，普通台式机和笔记本都能运行
- 输出格式多样

我就看中它一点，速度快。

软件安装异常的简单，因为提供了预编译的64位可执行文件

```bash
wget http://github.com/bbuchfink/diamond/releases/download/v0.9.25/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
# 有root全新啊
sudo mv diamond /usr/local/bin
# 无root权限, ~/bin是自己当前目录下
mv diamond ~/bin
```

因为 diamon的功能就是将蛋白或者翻译后的核苷酸和蛋白数据库进行比对，没有BLAST那么多功能，所以软件使用也是异常的简单。

第一步: 先从NCBI上下载蛋白数据库。 NR库是NCBI的非冗余蛋白数据库，

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz
```

也可以从[ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/](ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/)下载植物的蛋白数据库

第二步: 建库。就两个参数，`--in nr`输入文件，`--db nr` 输出的数据库前缀. 氨基酸序列中的结尾可以有"*"

```bash
diamond makedb --in nr --db nr
```

**注**: 假如要根据GFF提取蛋白序列，一定要注意输出的氨基酸序列中不能有"."在序列中，否则会报错。可以通过`seqkit grep -s -vrp '"\."' input.fa > output.fa` 进行过滤。

第三步: 搜索。就两个子命令，blastp和blastx，前者比对蛋白，后者比对DNA序列

```bash
diamond blastx --db nr -q reads.fna -o dna_matches_fmt6.txt
diamond blastp --db nr -q reads.faa -o protein_matches_fmt6.txt
```

`-q/--query`输入检索序列，`--out/-o`输出文件，默认以`--outfmt 6`输出结果和BLAST+的`--outfmt 6`结果一致。

注意事项:

- 默认参数主要是针对短序列，对于比较长的序列，使用`--sensitive`或`--more-senstive`提高敏感度。
- 默认的e-value阈值是0.001, 而BLAST是10，因此会比BLAST结果更加严格

性能优化:

- 设置比较低的`-e`参数
- 设置`-k`参数，减少输出的联配数目。这会降低临时文件大小和最终结果
- `--top`会输出得分比最好的分数低一定百分比的结果，
- `--compress 1`: 输出结果会以gzip进行压缩

## 参考文献

Benjamin Buchfink, Chao Xie, and Daniel H. Huson. Fast and sensitive protein alignment
using diamond. Nature methods, 12(1):59–60, Jan 2015.