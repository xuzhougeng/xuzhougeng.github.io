---
title: 如何安装JCVI
date: 2020-10-22 02:23:07.588
updated: 2020-11-02 19:21:23.927
url: /archives/how-to-install-jcvi
categories: 生信软件工具箱
tags: 软件安装
---

使用conda安装

```Bash
conda create -y -c bioconda -n jcvi jcvi 
```

使用时需要启动环境

```Bash
conda activate jcvi
```

在此基础上可以更新到最新版的JCVI (Github版本会不断的增加新功能并修订bug）

```Bash
pip install git+git://github.com/tanghaibao/jcvi.git
```

安装额外的依赖环境

- Kent tools

- BEDTOOLS

- EMBOSS

- [LAST](http://last.cbrc.jp/)

- LaTex

依赖环境中，BEDTOOLS, EMBOSS和LAST可以用conda安装

```Bash
# BEDTOOLS 和 EMBOSS
conda install -y -n jcvi -c bioconda bedtools emboss last

```

Kent tools需要编译, 并且依赖MySQL，暂时不需要安装。

```Bash
# Kent tools, 134M
wget http://hgdownload.cse.ucsc.edu/admin/jksrc.zip 
unzip  jksrc.zip
cd kent/src/lib
make 
```

对于Latex，conda的安装版本 ( texlive-core )在我后续的测试中各种出问题, 我最后根据根据官方的[quickinstall.html](https://www.tug.org/texlive/quickinstall.html) 进行了安装，具体的安装过程见 [LaTex](/archives/install-latex-without-root)

建议通过root进行安装（或者手工安装）。

```Bash
# ubuntu
sudo apt-get install -y texlive texlive-latex-extra texlive-latex-recommended
# centos
sudo yum install -y  texlive texlive-latex texlive-xetex texlive-collection-latexrecommended
```

最终通过Python版的MCscan对安装结果进行测试 （需要用conda安装seqkit)

```Bash
mkdir -p test && cd test
#Athaliana
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-44/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-44/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-44/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.44.gff3.gz
#Osativa
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-44/fasta/oryza_sativa/cdna/Oryza_sativa.IRGSP-1.0.cdna.all.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-44/fasta/oryza_sativa/pep/Oryza_sativa.IRGSP-1.0.pep.all.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-44/gff3/oryza_sativa/Oryza_sativa.IRGSP-1.0.44.gff3.gz 
# convert gff to bed
python -m jcvi.formats.gff bed --type=mRNA --key=transcript_id Arabidopsis_thaliana.TAIR10.44.gff3.gz > ath.bed
python -m jcvi.formats.gff bed --type=mRNA --key=transcript_id Oryza_sativa.IRGSP-1.0.44.gff3.gz > osa.bed 
# deduplication
python -m jcvi.formats.bed uniq ath.bed
python -m jcvi.formats.bed uniq osa.bed 
# Athaliana
seqkit grep -f <(cut -f 4 ath.uniq.bed ) Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz | seqkit seq -i > ath.cds
seqkit grep -f <(cut -f 4 ath.uniq.bed ) Arabidopsis_thaliana.TAIR10.pep.all.fa.gz | seqkit seq -i > ath.pep 
# Osativa
seqkit grep -f <(cut -f 4 osa.uniq.bed )  Oryza_sativa.IRGSP-1.0.cdna.all.fa.gz | seqkit seq -i  > osa.cds
seqkit grep -f <(cut -f 4 osa.uniq.bed ) Oryza_sativa.IRGSP-1.0.pep.all.fa.gz | seqkit seq -i  > osa.pep 
# create a new directory
mkdir -p cds && cd cds
ln -s ../ath.cds ath.cds
ln -s ../ath.uniq.bed ath.bed
ln -s ../osa.cds osa.cds
ln -s ../osa.uniq.bed osa.bed 
# compara
python -m jcvi.compara.catalog ortholog --no_strip_names ath osa 
```

如果顺利，当前目录下会有一个 ath.osa.pdf 文件

我测试的时候遇到了两个不顺利，都是和LaTex有关的报错。

第一个是latex没有安装导致的报错

```Bash
FileNotFoundError: [Errno 2] No such file or directory: 'latex': 'latex'
RuntimeError: Failed to process string with tex because latex could not be found

```

第二个是用conda安装texlive-core后，出现的第二个问题，还是得手工安装LaTex才解决。

```Bash
kpathsea: Running mktexfmt latex.fmt
Can't locate mktexlsr.pl in @INC  
```

