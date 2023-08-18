---
title: 从零开始学Circos绘制圈图(一)
date: 2019-07-26 13:45:44.293
updated: 2019-09-02 10:37:00.419
url: /archives/从零开始学Circos绘制圈图(一)
categories: 生信软件工具箱
tags: 可视化 | CIRCOS
---

一般基因组文章都会有下面这种酷炫图，用来描述基因组的基因密度分布，转座子的密度分布，和其他物种或者多倍体的多套染色体间的共线性关系，以及其他各种你只要测序就能加上的信息，比如说你要是测了ATAC-seq，加上全基因组开放状态，要是测了多个组织，多个时期的RNA-seq，那就加上热图展现这种变化关系。

![circos绘制基因组](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/7/1564118925560-3386fb85ab5a48d8b92b9635a30563e4.png)

当然除了基因组文章，其他类型的文章也可以考虑这种图。接下来我将会写一些列教程（可能有视频），通过教别人学Circos的方式来自学Circos。

## 环境配置

建议在Linux环境下配置Circos，之后只要用conda就能配置好分析环境

```bash
# 安装
## circos
conda create -c bioconda -n circos circos
```

测试软件安装结果

```bash
# 测试circos
conda activate circos
# 确认安装
circos -V
# 显示如下
# circos | v 0.69-8 | 15 Jun 2019 | Perl 5.026002
```

可以从<http://circos.ca/software/download/>下载官方的教程文件，分别是

- [circos-course-2017.tgz](http://circos.ca/distribution/circos-course-2017.tgz) 
- [circos-tutorials-current.tgz](http://circos.ca/distribution/circos-tutorials-current.tgz)

## 处理过程

Circos依赖于一些列的配置文件，用来定义复杂图形的各个部分，最终加工成图形。

因此，用Circos画图是一个不断增添内容的过程，你要不断根据输出结果来调整输入参数。

并且整个分析中，你还要拥有过关的数据预处理的能力，这是因为Circos不是数据处理工具，它只是展示你已有的数据。

![处理过程](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/7/1564114612076-b1639241644b45c8b802a39659340ed9.png)

## 快速开始

不管怎么样，先快速绘制出一个Circos图再说。

> 果子老师说过，我们不是先成为了老司机才开车，而是开车多了才成为了老司机。

第一步，先新建一个文件夹，用于存放本次分析的所有数据和配置文件

```bash
mkdir -p my_first_circos && cd my_first_circos
```

然后用`vim karyotype.tair10.txt`编辑文本，新增如下内容

```bash
chr - chr1 chr1 0 30427617 black
chr - chr2 chr2 0 19698289 black
chr - chr3 chr3 0 23459830 black
chr - chr4 chr4 0 18585056 black
chr - chr5 chr5 0 26975502 black
```

之后创建一个`circos.conf`文件，用于增加各类配置参数

```bash
touch circos.conf
```

用`vim circos.conf`，增加我们的第一条记录，染色体信息

```bash
karyotype = karyotype.tair10.txt
```

然而要想真正的出图，还需要增加至少以下配置语句才行

```bash
<ideogram>
<spacing>
default = 0.005r
</spacing>
radius           = 0.90r
thickness        = 20p
fill             = yes
stroke_color     = dgrey
stroke_thickness = 2p
</ideogram>

<image>
<<include etc/image.conf>>
</image>

<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
```

在当前路径下运行`circos -conf circos.conf`, 最终效果图如下

![第一张图](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/7/1564116740621-63c9fc5b9c0c4cddafd3eb7763653642.png)

虽然图比较丑，但是至少我们成功运行了人生第一次的circos, 这就相当于买了一套毛坯房，后面要做的事情就是不断装修。

比如说，我们至少可以让不同染色体拥有不同的颜色，修改之前的`karyotype.tair10.txt`中的最后一列

```bash
chr - chr1 chr1 0 30427617 chr1
chr - chr2 chr2 0 19698289 chr2
chr - chr3 chr3 0 23459830 chr3
chr - chr4 chr4 0 18585056 chr4
chr - chr5 chr5 0 26975502 chr5
```

在当前路径下运行`circos -conf circos.conf`, 最终效果图如下

![1564117790558](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/7/1564117790558-25d1cec4a31e4dd1beecde09134b1c98.png)

这就引出了第一个知识点，**配色**。

为了实现配色，需要`circos.conf`文件了有一个和配色有关的语句

```bash
<<include etc/colors_fonts_patterns.conf>>
```

这里`<<>>`表示通过**相对路径**的方式加载另外一个配置文件，它的实际路径是和`circos`所在目录同级的`etc`,可用下面语句看到`colors_fonts_patterns.conf`的内容

```bash
circos_path=$(dirname `which circos`)
less ${circos_path%bin}/etc/colors_fonts_patterns.conf
```

你会发现，这个文件里还嵌套其他的配置文件。最终通过层层排查，你才知道`etc/colors.ucsc.conf`才是实际定义我们填写的颜色名的文件，而颜色的定义如下：

```bash
chr1  = 153,102,0
chr2  = 102,102,0
chr3  = 153,153,30
chr4  = 204,0,0
chr5  = 255,0,0
```

还有一个问题，为什么这里用的是两个尖括号`<<`，而不是一个尖括号`<`呢？这是因为`<`已经被用于分隔不同的语句块，如下语句就表示`etc/image.conf`里的配置信息是用来调整和`image`有关的配置，而不是去调整`ideogram`的配置。

```bash
<image>
<<include etc/image.conf>>
</image>
```

---

以上是快速开始部分，后续将会在此基础上，做出发表级别的图。