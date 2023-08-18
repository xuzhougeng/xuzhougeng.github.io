---
title: 使用LACHESIS利用HiC数据进行基因组组装
date: 2019-11-01 23:44:28.716
updated: 2019-11-01 23:44:28.716
url: /archives/scaffolding-contig-with-LACHESIS
categories: 生信软件工具箱
tags: 组装 | Hi-C
---


![LACHESIS](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/11/image-540cd58f9f79426ea648f39202d1b2d5.png)

LACHESIS是一个比较早利用HiC数据辅助基因组组装的工具，文章发表在[Nature Biotechnology](https://doi.org/10.1038/nbt.2727)上.但是这款软件在2年前就不在维护了，GitHub主页上也推荐使用<https://github.com/theaidenlab>里的工具进行组装。

不过目前依旧有很多公司在用这个软件做基因组组装，所以我也想试试看。

## 软件安装

LACHESIS需要安装在Linux系统系统，最低内存不低于16GB，要求最小堆栈大小为10MB，可以用`ulimit -s`查看，通过`ulimit -s 10240`。同时依赖于以下软件

- gcc
- zlib
- boost C++ : > 1.52.0, < 1.67.0, 老版本只有1.57能被我下载
- SAMtools: < 0.1.19

以没有root权限为例，安装SAMtools和Boost这两个依赖库

```bash
mkdir ~/src
mkdir -p ~/opt/{biosoft,sysoft}
# 安装boost C++, 这一步会很久
cd ~/src
wget http://sourceforge.net/projects/boost/files/boost/1.57.0/boost_1_57_0.tar.bz2
tar xf boost_1_57_0.tar.bz2
cd boost_1_57_0
./bootstrap.sh --prefix=~/opt/sysoft/boost_1_57_0 --with-libraries=all --with-toolset=gcc
./b2 toolset=gcc
./b2 install
./bjam install
# 安装SAMtools
cd ~/src
wget https://astuteinternet.dl.sourceforge.net/project/samtools/samtools/0.1.19/samtools-0.1.19.tar.bz2
tar xf samtools-0.1.19.tar.bz2
cd samtools-0.1.19
make -j 20
```

然后是下载并解压缩

```bash
cd ~/opt/biosoft/
curl -o LACHESIS.zip https://codeload.github.com/shendurelab/LACHESIS/legacy.zip/master
unzip LACHESIS.zip
mv shendurelab-LACHESIS-2e27abb LACHESIS
cd LACHESIS
```

由于之前安装的samtools并不是安装在`/usr`目录下，因此要修改`src/include/gtools/`目录下的`SAMStepper.h`和`SAMStepper.cc`中的`#include <bam/sam.h>`, 指向实际地址。例如我的是`#include "/home/xzg/src/samtools-0.1.19/sam.h"`

```bash
export LACHESIS_BOOST_DIR=$HOME/src/boost_1_57_0
export LACHESIS_SAMTOOLS_DIR=$HOME/src/samtools-0.1.19
./configure --with-samtools=$HOME/src/samtools-0.1.19 \
    --with-boost=$HOME/src/boost_1_57_0/ 
    
./configure --with-samtools=/data/src/samtools-0.1.19 \
    --with-boost=/data/src/boost_1_57_0/ 
make -j 20
mv src/bin .
mv src/Lachesis bin/
```

> 如果服务器管理员帮你安装了一个高版本的boost，你需要让他进行卸载，否则默认会线用系统的boost，然后就会出错。关于更多的报错，参考[HiC软件安装篇之Lachesis](https://www.jianshu.com/p/51dc11d57ab4)

如果你Perl版本不是 5.14.2 ，那么我们需要打开bin下面的perl脚本，删除如下信息

```bash
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This software and its documentation are copyright (c) 2014-2015 by Joshua //
// N. Burton and the University of Washington.  All rights are reserved.     //
//                                                                           //
// THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS  //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF                //
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  //
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY      //
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT //
// OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR  //
// THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
```

而且还得将第一行替换成`#!/usr/bin/perl -w`

之后需要将LACHESIS加入环境变量

```bash
export PATH=$PATH:~/biosoft/LACHESIS/bin 
```

为了保证软件在后续不会出错，我们还要测试下

```bash
cd bin
./Lachesis INIs/test_case.ini
```

当你看到最后结果是`: Done!`就表明运行成功了，最终结果输出在当前目录下`outs`中

## 软件使用

实际的软件使用要求你需要提供至少两类输入文件

- HiC数据，双端fastq格式
- 初始组装，fasta格式

参考[使用ALLHiC基于HiC数据辅助基因组组装](/archives/Assemble-genome-using-ALLHiC-with-HiC-Data)的第一步，第二步和第三步，使用bwa-aln + bwa-sampe将fastq回帖到参考基因组，然后对数据进行过滤。那么在后续的流程都假设你有如下两个文件

- draft.asm.fasta: 初步组装结果
- sample.clean.bam: HiC数据比对的预处理结果

然后我们需要拷贝一个配置运行文件

```bash
cp ~/opt/biosoft/LACHESIS/bin/INIs/test_case.ini sample.ini
```

接着去修改参数，由于每个参数都会有说明，这里就展示我编辑过的几个参数

```bash
SPECIES = test # 写实际的物种名即可
DRAFT_ASSEMBLY_FASTA = draft.asm.fasta # 待组装序列的实际位置
SAM_DIR = . #表示当前目录下查找文件
SAM_FILES = sample.clean.bam #bam文件名
RE_SITE_SEQ = AAGCTT #酶切识别序列
USE_REFERENCE = 0 #不使用参考序列
BLAST_FILE_HEAD = . # BLAST的输出结果
CLUSTER_N = 16 # 最终聚类数目
```

上面这些参数的修改属于基本操作，而下面的参数可能需要反复修改

```bash
# contig中最小的酶切位点数,
CLUSTER_MIN_RE_SITES=25
# contig中最大的link密度, 也就是一个区域与多个contig存在信号
# 可能是异染色质或重复序列
CLUSTER_MAX_LINK_DENSITY=2
# 对于CLUSTER_MIN_RE_SITES过滤掉的contig在初步聚类后还有机会加入已有的分组中
# 如果它加入其中的信号是加入另一组信号的3倍
CLUSTER_NONINFORMATIVE_RATIO=3
# 允许成组的最小酶切数
ORDER_MIN_N_RES_IN_TRUNK=15
```

最终运行即可

```bash
ulimit -s 10240
Lachesis sample.ini
```

我们来构建最终的组装fasta

```bash
CreateScaffoldedFasta.pl draft.asm.fasta out/test_case
```

最终结果输出在`out/test_case/Lachesis_assembly.fasta`

使用体验:

- 软件安装有很多问题
- 无法处理多倍体
- 已经许久不维护
- 尝试了自己的物种，结果在排序步骤出现了问题

> 目前的想法: 这个软件我不会再用第二次了。