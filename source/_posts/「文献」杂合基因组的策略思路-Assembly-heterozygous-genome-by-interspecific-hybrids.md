---
title: 「文献」杂合基因组的策略思路
date: 2019-10-16 18:48:45.776
updated: 2019-10-24 10:13:03.463
url: /archives/Assembly-heterozygous-genome-by-interspecific-hybrids
categories: 文献阅读
tags: 组装
---

文献出处: [Sequencing a Juglans regia × J. microcarpa hybrid yields high-quality genome assemblies of parental specie](https://www.nature.com/articles/s41438-019-0139-1)

文章的亮点在于通过对一个F1子代进行三代测序，之后利用BioNano组装出两个亲本的光学图谱，最后根据光学图谱从F1中拆分出两套单倍型。

杂合基因组组装的时候，通常会通过构建近交系，花粉加倍(DH)来降低杂合度。文章采用了种间杂交的方式来避免杂合。这是因为种间杂交的基因组通常包括亲本的两个单倍型，因此能够从中分离出两个亲本的单倍型基因组。种间杂交的一个优势在于，它比较容易构建。

分析策略流程图如下

![流程图](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/10/1571214701466-aa197a4f932a4d26b1283283d7604c63.png)
文字版讲解如下：

Step1: 使用两个近源物种（都有基因组草图）进行杂交得到F1(MS1-56)，使用PacBio进行测序组装（460条contig, 1.05G, 两倍基因组大小）

Step2: 使用BioNano构建F1和两个亲本的光学图谱。

Step3: 通过自我比对（self-alignment），鉴定haploid和diploid。 亲本的光学图谱包含两个部分，一个是haploid（两个单倍体因为相似而坍缩成一个)，另一个是diploid(两个单倍体因为不够相似被拆分， 也就是分型(phased))。大概是下面这个效果，a表示Serr ctg4上存在两套单倍型的情况.

![光学图谱自我比对](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/10/1571221028944-3d69ebf2957849ab8a0c8fb1955bd2a5.png)

Step4: 将contig回帖到杂合的光学图谱上，组装成scaffolds，大小为1.06G。之后根据杂合后代光学图谱和亲本光学图谱的关系，对scaffold进行拆分。上图b,c表示的是将F1后代的光学图谱和双亲的光学图谱进行比较，来鉴定出子代中属于一方亲本的单倍型。

Step5: 上述因为太短而无法回帖到光学图谱的contig，比对到illumina组装的双亲基因组，进行区分。

Step6: 根据高密度遗传图谱对scaffold进行排序和确认方向，最终得到准染色体级别的基因组。

最后基因组还可以补洞和纠错，来提高质量，比如说10X Genomics Chicago 测序。

> 最后一点感想：文章利用光学图谱进行分型的思路的确很棒。但是我想的是，为啥不直接测亲本基因组，用HiC进行分型呢？