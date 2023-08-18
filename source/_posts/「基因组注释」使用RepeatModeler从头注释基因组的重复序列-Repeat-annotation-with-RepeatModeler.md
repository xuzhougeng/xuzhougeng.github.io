---
title: 「基因组注释」使用RepeatModeler从头注释基因组的重复序列
date: 2019-08-26 17:31:51.898
updated: 2019-08-27 22:32:24.936
url: /archives/Repeat-annotation-with-RepeatModeler
categories: 生信软件工具箱
tags: 注释 | 重复序列
---

RepeatModeler可用来从头对基因组的重复序列家族进行建模注释，它的核心组件是RECON和RepatScout。

## 使用方法

以拟南芥的参考基因组为例，假设基因组的名字为"Athaliana.fa"

第一步：为RepeatModeler创建索引数据库

```bash
BuildDatabase -name ath -engine ncbi Athaliana.fa
# -engine ncbi： 表示使用rmblast
# -name aht： 表示数据库的名字为ath
```

这一步其实认为调用了`makeblastdb`，结果文章和`makeblastdb`一致

第二步：运行RepeatModeler

```bash
RepeatModeler -database ath -engine ncbi -pa 10 &> ath.out &
# -database 要和上一步一致
# -engine 要和上一步一致
# -pa 表示线程数
```

运行中的的文件存放在`RM_.xxx`文件夹下

```bash
RM_82213.MonOct151054032018
├── consensi.fa 
├── consensi.fa.classified
├── consensi.fa.masked
├── families-classified.stk
├── families.stk
├── round-1
├── round-2
├── round-3
├── round-4
└── round-5
```

运行结束后，就得到了`ath-families.fa`和`ath-families.stk`。 前者是找到的重复序列，后者是Stockholm格式的种子联配文件(seed alignment file), 可以用`util/dfamConsensusTool.pl`上传到`Dfam_consensus`数据库中。

```bash
grep -i Unkown ath-families.fa
# 没有结果，全都归类，毕竟拟南芥
```

`RM_.xxx`的`consensi.fa.classified`和`ath-families.fa`内容一样，也是FASTA格式的文件，，只不过每个序列的ID会标注它来自于哪个重复序列家族，如果无法归类，则用"Unkown"标注。

你可以将`ath-families.fa`分ModelerID.lib和Modelerunknown.lib, 

```bash
seqkit grep -nrp Unknown ath-families.fa > Modelerunknown.lib
seqkit grep -vnrp Unknown ath-families.fa >  ModelerID.lib
```

其中Modelerunknown.lib用RepeatMasker或TEclass进一步注释，如果能够被分类则从Modelerunknown.lib,移动到ModelerID.lib

这里用TEclass进行分类，软件安装使用参考[使用TEclass对TE一致性序列进行分类](/archives/Classify-TE-consensus-sequence-by-TEclass).

```bash
TEclassTest Modelerunknown.lib
```

> 这一步其实无法运行，因为没有输入。

第三步：过滤基因片段

> 许多文章都没有这样子做，因此这一步完全可以不用看

RepeatModeler找到的重复序列进一步在植物蛋白数据库（不包括转座子蛋白）进行搜索，如果和植物蛋白匹配，或者在序列的侧翼50bp以内，就将该重复序列剔除。这一步可用工具是`ProtExcluder`
ProtExcluder的安装需要HMMER3， 方法如下

首先一定要装HMMER3.0，然后安装方法为

```bash
wget http://eddylab.org/software/hmmer/hmmer-3.0.tar.gz
tar xf hmmer-3.0.tar.gz
cd hmmer-3.0
./configure --prefix=$HOME/opt/biosoft/hmmer-3.0
make && make install
cd easel
./configure --prefix=$HOME/opt/biosoft/hmmer-3.0
make && make install
```

后续安装的ProtExcluder1.1非常坑爹，他居然认为hmmer的运行文件是放在binaries目录下，所以你还需要去`~/opt/biosoft/hmmer-3.0`把bin文件夹改成binaries

然后安装ProtExcluder1.1

```bash
tar zxf  ProtExcluder1.1.tar.gz
mv ProtExcluder1.1 ~/opt/biosoft
cd  ~/opt/biosoft/ProtExcluder1.1 
./Installer.pl   -m   ~/opt/biosoft/hmmer-3.0 -p P ~/opt/biosoft/ProtExcluder1.1/
# 注意运行路径后一定要有"/"
```

然后将找到的重复序列用blastx比对到植物的蛋白数据库中，你可以到RefSeq上进行下载。

```bash
blastx -query ath-families.fa -db /path/to/refseq-plant/db/plant.protein -num_threads 70 > ath-families.blastx &
```

最后运行ProtExcluder.pl

```bash
/path/to/ProtExcluder.pl -f 50   ath-families.blastx  ath-families.fa 
```

从理论上说结果是"XXXnoProtFinal"，但这个破软件各种报错，而且报错信息信息量太少，直接放弃这个破软件了。

考虑到目前看了很多文献也没人说要过滤，所以这一步就放弃吧

注1： 一般而言RepeatModeler在一周内就能运行完，物种小一点，服务器好一点，基本上一天就行了
注2： 一般而言，同源注释所用的重复序列库的数据比较可靠，但不一定全，而我们自己建立的重复序列库比较全，但是未必准。
注3： SwissProt的植物蛋白数据库的下载地址为：<ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions>，而NCBI的Refseq植物蛋白数据库下载地址为<ftp://ftp.ncbi.nih.gov/refseq/release/plant/>
