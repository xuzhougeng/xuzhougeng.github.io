---
title: 如何绘制物理图谱和遗传图谱的对应关系
date: 2021-12-13 06:06:48.219
updated: 2021-12-13 06:06:48.219
url: /archives/draw-match-between-physical-and-genetic-map
categories: 基因组学
tags: JCVI
---

唐海宝老师开发的JCVI有一个工具，叫做ALLMAPS， 能够展示遗传图谱和物理图谱的对应关系，如下所示

![alignment](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/12/image-095b4255703740339a37cb6b768e356e.png)

但是这个图的目标是为了对ALLMAPS的scaffold结果进行可视化，并不是专门用于展示遗传图谱的标记和物理图谱的对应关系。尽管在allmaps这个组件下提供了plot函数，命令行输入只要求 input.bed 和 seqid， 但实际运行的时候还要求 allmaps path的中间文件, xxxx.lifted.bed, xxxx.agp, weight.txt等文件。

为了解决这一问题，我阅读了allmaps.py的源代码，在plot的基础上增加了一个plot2函数，只需要用户输入 input.bed, 染色体编号和染色体的长度就能够画图。

```bash
# python -m jcvi.assembly.allmaps plot2 input.bed 染色体编号 染色体长度
python -m jcvi.assembly.allmaps plot2 input.bed chr1 123140023
```

其中input.bed的格式要求有6列，

- 标记所在染色体名
- 标记所在染色体的start
- 标记所在染色体的end, 通常就是start+1
- 标记对应的图谱位置， 要求输入为"图谱名-连锁图谱所在组:连锁图谱的遗传距离"
- 标记名

案例

```text
Chr1    68185909        68185910        male-14:48.470000        Chr1:68185910
Chr1    68479621        68479622        male-14:49.380000        Chr1:68479622
Chr1    68595299        68595300        male-14:48.440000        Chr1:68595300

```

可以先生成如下的csv文件，然后转换成bed

![input.csv](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/12/image-6dcf0d4e3efd45de994bfc1a7f646e97.png)

```python
python -m jcvi.assembly.allmaps merge  male.csv male.bed
```

目前有plot2函数的代码还在我的项目下，[xuzhougeng/jcvi](https://github.com/xuzhougeng/jcvi), 待代码稳定了，再PR。