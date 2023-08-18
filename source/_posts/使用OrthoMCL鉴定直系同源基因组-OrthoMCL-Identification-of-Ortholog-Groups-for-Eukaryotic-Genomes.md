---
title: 使用OrthoMCL鉴定直系同源基因组
date: 2019-09-03 15:57:48.135
updated: 2019-09-04 15:26:58.703
url: /archives/OrthoMCL-Identification-of-Ortholog-Groups-for-Eukaryotic-Genomes
categories: 生信软件工具箱
tags: 基因家族
---



OrthoMCL是目前最常用的基因家族分析软件，从2013年发布2.0版本之后再也没有更新过，虽然它的安装过程复杂负责，但是依旧挡不住大家对他的喜好。

当然软件安装复杂是相对于之前，现在用Docker就可以轻松的安装和使用OrthoMCL。由于

```bash
docker pull jasonkwan/orthomcl_docker
```

由于该容器自带一个`run_orthomcl.py`, 只要你准备好输入数据和配置文件，就能够自动化进行分析.

其中run_orthomcl.py的参数有5个

- -t/--table_path fasta_table的文件路径
- -s/--start_stage: 开始阶段，1,2,3,4
- -e/--end_stage: 结束阶段， 1,2,3,4
- -p/--processors: 线程数
- -c/--config_file: 配置文件

其中-s和-e的1-4分别对应为

- 1: 使用orthomclAdjustFasta, orthomclFilterFasta进行数据预处理
- 2: 使用DIAMOND进行all-vs-all blast
- 3: 使用MySQL数据库寻找配对
- 4: 使用MCL处理配对

## 分析流程

我们以一个具体的案例说明下。新建一个文件夹，例如说`test`， 里面是收集好的物种氨基酸序列。

```bash
$ ls
Alyrata.faa  Athaliana.faa  Chirsuta.faa
```

由于Docker版本的OrthoMCL使用DIAMOND进行blastp比对，因此一定要保证你的氨基酸序列中没有"."，否则会报错。可以用`seqkit grep`过滤不合格的序列。

```bash
seqkit grep -s -vrp '"\."' input.fa > output.faa
```

之后创建一个fasta_table，放在`test`目录下。该文件分为两列，第一列是文件名，第二列是缩写。

```bash
Alyrata.faa	Aly
Athaliana.faa	Ath
Chirsuta.faa	Chi
```

准备配置文件`orthomcl.config`，同样放在`test`目录下

```bash
dbVendor=mysql
dbConnectString=dbi:mysql:orthomcl:mysql_local_infile=1
dbLogin=root
dbPassword=PAssw0rd
similarSequencesTable=SimilarSequences
orthologTable=Ortholog
inParalogTable=InParalog
coOrthologTable=CoOrtholog
interTaxonMatchView=InterTaxonMatch
percentMatchCutoff=50
evalueExponentCutoff=-5
oracleIndexTblSpc=NONE
```

在test同级目录下运行如下命令

```bash
docker run --privileged=true  --rm -it --volume $PWD/test:/outdir:rw jasonkwan/orthomcl_docker:latest bash
```

就会进入Docker的交互命令行，运行`run_orthomcl.py`

```bash
cd outdir
run_orthomcl.py --table_path fasta_table --config_file orthomcl.config --processors 32 &
```

假如不希望进入交互命令行，那么需要按照下面的方法进行运行

```bash
docker run --privileged=true  --rm --volume $PWD/test:/outdir:rw jasonkwan/orthomcl_docker:latest /bin/bash -c "/tmp/.runconfig.sh && run_orthomcl.py --table_path /outdir/fasta_table --config_file /outdir/orthomcl.config --processors 32 "
```

最终输出结果是groups.txt。下一个问题是，当你有了groups.txt后，下面能做什么分析呢？

按照对基因组学文章的整理，基本上就是下面两个

- 使用单拷贝基因构建系统发育树。
- 使用CAFE进行基因家族扩张收缩分析

我正在找一篇文章尝试重现这个流程。

## 参考资料

- https://hub.docker.com/r/jasonkwan/orthomcl_docker
- https://bitbucket.org/jason_c_kwan/orthomcl_docker/src/master/Dockerfile
- 特别感谢[Leo](https://github.com/leoatchina)在Docker上的帮助
- Dockerfile中CMD和ENTRYPOINT: https://www.cnblogs.com/sparkdev/p/8461576.html
