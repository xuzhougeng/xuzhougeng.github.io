---
title: 使用TEclass对TE一致性序列进行分类
date: 2019-08-27 21:56:48.413
updated: 2022-07-19 07:52:08.69
url: /archives/classify-te-consensus-sequence-by-teclass
categories: 生信软件工具箱
tags: 注释 | 重复序列
---

## 软件安装

软件地址在<http://www.compgen.uni-muenster.de/tools/teclass/index.hbi?>, 由于TEclass这个软件已经许久没有更新了，因此还要讲解下安装步骤。

> 最近更新了一个2.1.3c, 经过我测试发现，应该就是把之前无法下载URL做了更新。classifiers.tar.gz无更新。

```bash
wget http://www.compgen.uni-muenster.de/tools/teclass/download/TEclass-2.1.3.tar.gz
tar xf TEclass-2.1.3.tar.gz
cd TEclass-2.1.3
```

下载依赖的软件

```bash
sh Download_dependencies.sh
```

 由于代码老旧，部分内容无法自动下载，需要手动下载， 例如librf, blast.  最终要保证文件夹下有如下文件

- libsvm.tar.gz: <http://www.csie.ntu.edu.tw/~cjlin/libsvm/>
- glimmer.tar.gz:  <http://ccb.jhu.edu/software/glimmer/>
- librf.tar.gz: <http://mtv.ece.ucsb.edu/benlee/librf.html>
- lvq_pak.tar: <http://www.cis.hut.fi/research/som-research/nnrc-programs.shtml>
- blast.tar.gz: <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED>

例如blast

```bash
curl -o 'blast.tar.gz' ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz
```

编译依赖的软件

```bash
sh Compile_dependencies.sh
```

![安装过程](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/8/2013053-3885005894b4ac70-0f20bdb0eaad4826a9dd4fdf81c587bc.png)

安装TEclass, 这一步可以跳过 RepBase的配置。

```bash
perl Configure.pl
```

安装预编译的分类器，这一步在TEclass的目录下运行，并解压缩

```bash
wget http://www.compgen.uni-muenster.de/tools/teclass/download/classifiers.tar.gz
mv classifiers.tar.gz classifiers/
cd classifiers
tar xf classifiers.tar.gz
```

测试运行

```bash
./TEclassTest.pl ./testfile.fa
```

## 软件使用

### 构建分类器

如果你想使用最新的RepBase，那么就需要自己从[RepBase](http://www.girinst.org/repbase/index.html)下载数据进行构建。

如果是单核处理器，可能要花费数周的时间，所以建议用上你的所有线程。

```bash
/TEclassBuild.pl -x 0  -o new_classifiers -p 99
```

### 重复序列分类

在RepeatModeler建模后，提取Unknown序列使用`TEclassTest`进行归类，假如输入文件命名为Modelerunknown.lib

```bash
TEclassTest Modelerunknown.lib
```

结果在`Modelerunknown.lib_xxxx`, `xxxx`是你运行日期。

```bash
Modelerunknown.lib # 输入文件
Modelerunknown.lib.html 
Modelerunknown.lib.lib # 输出结果
Modelerunknown.lib.stat #结果统计
``` 

`Modelerunknown.lib.lib`中的fasta会有分类信息，如

```bash
>rnd-1_family-12#Unknown ( RepeatScout Family Size = 705, Final Multiple Alignment Size = 88, Localized to 114 out of 117 contigs )|TEclass result: LTR|forward|ORFs: 583..2355:+1
```

> 需要注意的是，TEclass的输出结果是不被RepeatMasker所识别的，需要你更改原来的#Unknown为你的预测结果才行。

其他参数:

- `-c`: 训练的分类器所在路径, 默认是`TEclass-2.1classifiers`
- `-o`: 输出结果路径，默认在当前文件下新建
- `-r`: 预测输入序列的反向互补序列

参考文献:  TEclass: a tool for automated classification of unknown eukaryotic transposable elements