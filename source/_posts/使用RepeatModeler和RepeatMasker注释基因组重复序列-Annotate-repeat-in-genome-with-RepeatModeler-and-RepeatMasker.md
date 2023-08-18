---
title: 使用RepeatModeler和RepeatMasker注释基因组重复序列
date: 2019-09-18 10:58:40.865
updated: 2019-09-18 10:58:40.865
url: /archives/Annotate-repeat-in-genome-with-RepeatModeler-and-RepeatMasker
categories: 生信软件工具箱
tags: 重复序列
---

重复序列注释有两种常用策略，基于同源序列相似性和基于重复序列结构特征。`RepeatMasker`是基于同源序列相似性注释序列的常用工具, `RepeatModeler`可用来从头对基因组的重复序列家族进行建模注释，它的核心组件是RECON和RepatScout。

这篇教程介绍如何使用`RepeatModeler`从头鉴定基因组的重复序列，之后用`RepeatMasker`根据自定义的重复序列库注释基因组的重复序列。

## 软件安装

原本的`RepeatMasker`和`RepeatModeler`的手动安装需要配置很多文件，但是利用bioconda就只用一行命令。

```bash
conda create -n repeat repeatmasker repeatmodeler
```

之后在RepeatMasker环境下配置运行环境。由于我的miniconda装在`~/opt`路径下，因此对应的RepeatMasker路径为`~/opt/miniconda3/envs/repeat/share/RepeatMasker/`

```bash
conda activate repeat
cd ~/opt/miniconda3/envs/repeat/share/RepeatMasker/
perl ./configure
```

这一步只需要配置好比对软件

![配置比对工具](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/9/2013053-9aa560fcadf94885-cb7b3c03d6a644d5be1ebca74a1c1fa4.png)

之后就会显示RepeatMasker已经配置完毕，其中Dfam_3.0是用于注释的数据库。

![配置完成](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/9/1568772932341-2431026ea88645559eda6daf7ef96cd0.png)

## 软件运行

以拟南芥的参考基因组为例，基因组命名为"Athaliana.fa"

第一步：为RepeatModeler创建索引数据库

```bash
BuildDatabase -name ath -engine ncbi Athaliana.fa
# -engine ncbi： 表示使用rmblast
# -name aht： 表示数据库的名字为ath
```

第二步：运行RepeatModeler从头预测

```bash
RepeatModeler -database ath -engine ncbi -pa 20 &> ath.out &
# -database 要和上一步一致
# -engine 要和上一步一致
# -pa 表示线程数
```

这一步运行时间相对比较久，和线程数有关。运行中的的文件存放在`RM_.xxx`文件夹下

```bash
RM_100741.WedSep181006282019
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

`ath-families.fa`的fasta的序列部分格式为`>repeatname#class/subclass`，用于表明每个重复序列的归类。

第三步：根据自定义的重复序列数据库注释基因组

```bash
RepeatMasker -lib ath-families.fa -e ncbi -dir . Athaliana.fa
```

RepeatMasker比较常用的参数如下

- `-e`: 搜索引擎，默认都选择ncbi
- `-pa`: 并行计算，多线程
- `-s`, `-q`, `-qq`: 搜索速度，速度和敏感度成反比
- `-lib`: 自定义重复数据库
- `-species`: 指定物种，例如human, mouse, arabidopsis
- `-gff`: 额外输出GFF文件

输出结果中, 以.masked结尾的是用N屏蔽后的序列，以tal结尾的则是统计各种重复序列的比例。