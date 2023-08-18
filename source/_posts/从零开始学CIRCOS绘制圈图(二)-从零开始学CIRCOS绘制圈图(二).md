---
title: 从零开始学CIRCOS绘制圈图(二)
date: 2019-07-30 10:52:23.083
updated: 2019-09-02 12:59:45.903
url: /archives/从零开始学CIRCOS绘制圈图(二)
categories: 生信软件工具箱
tags: 可视化 | CIRCOS
---

在[从零开始学CIRCOS绘制圈图(一)](/archives/%E4%BB%8E%E9%9B%B6%E5%BC%80%E5%A7%8B%E5%AD%A6Circos%E7%BB%98%E5%88%B6%E5%9C%88%E5%9B%BE(%E4%B8%80))， 我们已经绘制出一个比较丑的circos图，这一部分是讲解一些细节。

这一部分会从上一步的两个文件开始，分别是

`karyotype.tair10.txt`

```bash
chr - chr1 chr1 0 30427617 chr1
chr - chr2 chr2 0 19698289 chr2
chr - chr3 chr3 0 23459830 chr3
chr - chr4 chr4 0 18585056 chr4
chr - chr5 chr5 0 26975502 chr5
```

`circos.conf`

```bash
karyotype = karyotype.tair10.txt

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

## 配色

第一个要解决的问题是配色问题，其实上一部分的配色定义依旧有点丑。如果不想重新打开文件修改的话，其实除了在karyotype文件里定义颜色外，我们还可以直接在circos.conf文件定义颜色.

举个例子，在`karyotype = karyotype.tair10.txt`后加一行

```bash
chromosomes_color = chr1=rdylbu-11-div-1,chr2=rdylbu-11-div-3,chr3=rdylbu-11-div-5,chr4=rdylbu-11-div-7,chr5=rdylbu-11-div-9
```

问题是，我们怎们知道这些名字后所代表的颜色呢？

Circos中颜色的命名格式为`PALETTE-NUMCOLORS-TYPE-IDX`:

- PALETTE:调色版名，如rdylbu
- NUMCOLORS: 颜色数目,  11
- 调色版类型: div(diverging), seq(sequential), qual(qualitative)
- IDX: 调色版中的颜色索引

而Circos颜色来自于[http://colorbrewer2.org](http://colorbrewer2.org/)

因此，`gnbu-9-seq`对应的是就是下图的`9-class GnBu`

![颜色](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/7/circos-colors-88bba063fa264c9f9825a9b9fd3baad0.png)

配色不仅仅用在karyotype中，在后续的热图和柱状图等还会涉及到它，毕竟一张好看的图，配色占了很大的比例。

## 显示标签

默认输出图片是没有染色体名字的标签，需要在`<ideogram>`里添加和`label`有关的参数

```bash
...
</spacing>
...
show_label     = yes #展示label
label_font     = default # 字体
label_radius   = dims(ideogram,radius) + 0.05r #位置
label_size     = 16 # 字体大小
label_parallel = yes # 是否平行

label_format   = eval(sprintf("%s",var(chr))) # 格式
</ideogram>
```

关于标签(label)， 更详细的介绍在<http://circos.ca/documentation/tutorials/ideograms/labels/>

重新运行之后，发现字相对而言太小了。有两种结局方案，一种是调整`label_size`，比如说48，另一种是调整图片整体大小。

![显示标签](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/7/circos-small-label-8fa0fec44a434954be54e07f8ce85ba8.png)

## 输出设置


默认情况下，输出图片的 半径是1500p, 所以设置的`label_size=16`就会显得特别小。我们可以在配置文件中调整图片的班级，以及其他设置。

```bash
<image>
dir*    = .    # 输出文件夹
radius* = 500p # 图片半径
svg*    = no   # 是否输出svg
<<include etc/image.conf>>
</image>
```

**注**: 在参数后加一个`*`表示覆盖已有的设置，比如说`svg*=no`就是覆盖已有的`svg=yes`。

运行结果如下

![修改输出图片大小](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/7/circos-large-label-3f6597855c3641b7baf918b7cc5fd558.png)

## 刻度(ticks)

大部分的circos还会有刻度来展示染色体的大小。由于ticks的定制比较复杂，所以一般会单独搞一个配置文件，`ticks.conf`存放ticks相关的参数设置.

```bash
chromosomes_units = 1000000
show_ticks          = yes
show_tick_labels    = yes

<ticks>

radius           = 1r
color            = black
thickness        = 2p

multiplier       = 1e-6 #输出的标签为实际长度与其相乘

format           = %d # %d表示显示整数

<tick>
spacing        = 1u
size           = 5p
</tick>

<tick>
spacing        = 5u
size           = 10p
show_label     = yes
label_size     = 10p
label_offset   = 10p
format         = %d
</tick>

</ticks>
```

用 `show_ticks`和 `show_tick_labels`控制是否展示刻度，以及刻度对应的标签。

`<ticks`与`</ticks>`里控制刻度的全局参数，例如位置为1r(radius=1r), 颜色为黑色(color=black), 厚度为2p(thickness=2p),   由于默认直接展示染色体的实际位置，因此会显示1000000这种结果，所以定义`multiplier=1e-6`,  实际显示结果为 1000000 * 1e-6 = 1。

后面就可以通过`<tick>`和 `</tick>`来分别绘制不同类型的tick，重要的参数如下：

- `spacing`表示刻度之间的距离，1u表示一个长度单位，需要在`circos.conf`文件里通过`chromosome_unit`来定义，通常都是`chromosome_unit=1000000`。
- `size`表示tick的长度
- `show_label`: 控制是否展示标签，默认不展示。
- `label_offset`: 则是让label往外在偏移一点距离

`circos -conf circos.conf`运行结果如下

![增加刻度](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/7/circos-ticks-c284af9cff484548983579dcee75c154.png)

我们发现有些刻度的标签和染色体标签发生了重叠，这个可以通过`label_radius`进行调整。

## 单位

上面出现了控制图形不同元素大小的三个单位，p,r,u。p(pixels), 表示绝对大小， r(relative), 相对大小， u(chromosome unit)。 如果使用p作为单位，需要考虑最终输出图形`<image>`定义的radius。 而r是相对大小，会随着最终图形大小而发生变换。u一般在显示刻度时使用。

---

这一部分是在原来简陋的输出上进行了美化，没有用到除了染色体长度以外的信息。下一部分介绍如何展示基因密度等信息。
