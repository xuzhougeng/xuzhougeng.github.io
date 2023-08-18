---
title: 使用refgenie管理你的参考基因组
date: 2019-09-09 20:41:45.089
updated: 2019-09-09 20:41:45.089
url: /archives/Manage-reference-genome-with-refgenie
categories: 生信软件工具箱
tags: 流程工具 | 服务器
---


在服务器管理初期，我管理参考基因组的方法非常简单，就是以物种名建立一个文件夹，然后把和该物种有关的FASTA文件，GFF文件都放到该文件中，之后在文件夹中建立不同软件的索引。

当然目前已经有一些软件可以帮助你进行管理，比如说`refgenie`。它是一个Python编写的参考基因组管理工具，软件的设计思路如下:

![设计思路](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/9/image-799c4ee8c367485080b144fe4bec4c74.png)

你既可以在本地自己构建，也可以从它的服务器上下载已有的物种。

## 软件安装

refgenie的安装非常简单，只需要一行代码

```bash
# 建议先安装miniconda
pip install refgenie
```

之后，你需要初始化一个文件夹，之后下载的基因组都会在该目录下

```bash
mkdir -p ~/reference/
refgenie init -c ~/reference/genome_config.yaml
```

将这一行`export REFGENIE=~/reference/genome_config.yaml`根据所用的SHELL加入到对应的`.bashrc`或`.zshrc`

## 软件使用

我们可以使用`refgenie listr`去查看目前refgenomes服务器中已经有的参考基因组，输出信息如下

```bash
Querying available assets from server: http://refgenomes.databio.org/assets
Remote genomes: hg19, hg19_cdna, hg38, hg38_cdna
Remote assets:
  hg19: bismark_bt1_index; bismark_bt2_index; bowtie2_index; bwa_index; fasta; hisat2_index
  hg19_cdna: bowtie2_index; hisat2_index; kallisto_index; salmon_index
  hg38: bismark_bt1_index; bismark_bt2_index; bowtie2_index; bwa_index; fasta; hisat2_index
  hg38_cdna: bowtie2_index; hisat2_index; kallisto_index; salmon_index
```

我们可以用`refgenie pull`来下载数据

```bash
refgenie pull --genome hg38 --asset bowtie2_index
```

当然更常见的情况是，你的物种并不在已有的列表中，以及这是一个国外服务器，你甚至都无法拉取列表，所以我们需要用`refgenie build`来构建参考基因组。

我们以拟南芥参考基因组为例，我们需要先从[EnsemblPlants](http://plants.ensembl.org/index.html)上下载参考基因组序列

```bash
wget 'ftp://ftp.ensemblgenomes.org/pub/plants/release-44/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz'
```

先导入fasta文件

```bash
refgenie build --genome TAIR10 --asset fasta --fasta Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
```

建立bwa的索引

```bash
refgenie build --genome TAIR10 --asset bwa_index --fasta Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
```

此时检查`~/reference/genome_config.yaml`文件，里面记录了刚才新建文件的位置

```bash
config_version: 0.2
genome_folder: /path/to/reference
genome_server: http://refgenomes.databio.org
genomes:
  TAIR10:
    assets:
      fasta:
        path: fasta/TAIR10.fa
        asset_description: Sequences in the FASTA format
      fai:
        path: fasta/TAIR10.fa.fai
        asset_description: Indexed fasta file, produced with samtools faidx
      chrom_sizes:
        path: fasta/TAIR10.chrom.sizes
        asset_description: Chromosome sizes file
      bwa_index:
        path: bwa_index
        asset_description: Genome index for Burrows-Wheeler Alignment Tool, produced with bwa index
```

和build相关的详细信息参考<http://refgenie.databio.org/en/latest/build/>

## 参考资料

官方文档: <http://refgenie.databio.org/en/latest/>![refgenie_interfaces.svg](1)