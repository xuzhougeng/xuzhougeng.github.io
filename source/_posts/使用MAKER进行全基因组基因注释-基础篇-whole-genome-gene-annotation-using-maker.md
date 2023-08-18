---
title: 使用MAKER进行全基因组基因注释-基础篇
date: 2019-08-28 09:22:20.738
updated: 2020-07-07 07:16:35.541
url: /archives/whole-genome-gene-annotation-using-maker
categories: 生信软件工具箱
tags: 注释 | MAKER
---

![maker](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/2013053-3219030babf92e35-42d10cc44f5d4d159cecec7bc04c0087.png)

在基因组注释上，MAKER算是一个很强大的分析流程。能够识别重复序列，将EST和蛋白序列比对到基因组，进行从头预测，并在最后整合这三个结果保证结果的可靠性。此外，MAKER还可以不断训练，最初的输出结果可以继续用作输入训练基因预测的算法，从而获取更高质量的基因模型。

Maker的使用比较简单，在软件安装成后，会有一个"data"文件夹存放测试数据

```bash
ls ~/opt/biosoft/maker/data
dpp_contig.fasta  dpp_est.fasta  dpp_protein.fasta  hsap_contig.fasta  hsap_est.fasta  hsap_protein.fasta  te_proteins.fasta
```

以"dpp"开头的数据集为例，protein表示是同源物种的蛋白序列，est是表达序列标签，存放的是片段化的cDNA序列，而contig则是需要被预测的基因组序列。

让我们新建一个文件夹，并将这些测试数据拷贝过来。

```bash
mkdir test01 ; cd test01
cp ~/opt/biosoft/maker/data/dpp* .
```

由于基因组注释设计到多个程序，多个步骤，每个步骤可能都有很多参数需要调整，因此就需要建立专门的配置文件用来告诉maker应该如何控制流程的运行。

如下步骤创建三个以ctl结尾的配置文件

```bash
~/opt/biosoft/maker/bin/maker -CTL
ls *.ctl
maker_bopts.ctl  maker_exe.ctl  maker_opts.ctl
```

- maker_exe.ctl: 执行程序的路径
- maker_bopt.ctl: BLAST和Exonerat的过滤参数
- maker_opt.ctl: 其他信息，例如输入基因组文件

maker\_exe.ctl和maker\_bopt.ctl可以简单用less查看，可不做修改，maker\_opt.ctl是主要调整的对象。 使用`vim maker_opt.ctl`修改如下内容

```bash
genome=dpp_contig.fasta
est=dpp_est.fasta
protein=dpp_protein.fasta
est2genome=1
```

修改完之后多花几分钟看看每个参数的设置，尽管很枯燥，但是考虑这个工具你可能会反复多次使用，所以这点时间是一定要花的。

随后就可以在当前路径运行程序

```bash
~/opt/biosoft/maker/bin/maker &> maker.log &
```

输出结果见"dpp_contig.maker.output", 重点是"dpp_contig_master_datastore_index.log"文件，由于maker会拆分数据集并行计算，因此该文件记录总体的运行情况，需要关注其中是否有"FAILED","RETRY","SKIPPED_SAMLL","DIED_SIPPED_PERMANET"，因为这意味着有些数据出于某些原因没有运算。

最后，我们需要将并行运算的结果进行整合，导出GFF文件, 转录本序列和蛋白序列

```bash
~/opt/biosoft/maker/bin/fasta_merge -d dpp_contig_master_datastore_index.log
~/opt/biosoft/maker/bin/gff3_merge -d dpp_contig_master_datastore_index.log
```

在该目录下就会出现, "dpp_contig.all.gff", "dpp_contig.all.maker.proteins.fasta","dpp_contig.all.maker.transcripts.fasta"

其中GFF文件就需要用IGV，JBrowse, Apollo下展示来检查下注释是否正确。


## 附录

**软件安装**：MAKER可以免费用于学术用途，但是未经许可不可商用。目前有两个版本2018年5月4日更新的2.31.10和测试版3.01.02.出于稳定性考虑，安装前者。后续假设已经在<http://yandell.topaz.genetics.utah.edu/cgi-bin/maker_license.cgi>进行登记，并且下载了压缩包"maker-2.31.10.tgz"

先检查下自己的系统情况，看需要补充哪些库

```bash
tar xf maker-2.31.10.tgz
cd maker/src
perl Build.PL
```

这一步之后会罗列出后续需要运行的命令来完成安装

```bash
./Build installdeps
./Build installexes
./Build install
./Build status
```

## 参考资料

- Genome Annotation and Curation Using MAKER and MAKER-P