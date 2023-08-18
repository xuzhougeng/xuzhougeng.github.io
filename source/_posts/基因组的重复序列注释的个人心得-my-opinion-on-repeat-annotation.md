---
title: 基因组的重复序列注释的个人心得
date: 2022-07-20 03:46:27.423
updated: 2022-07-21 02:12:56.625
url: /archives/my-opinion-on-repeat-annotation
categories: 基因组学
tags: 流程工具
---

> 2022-07-20 第一版
> 2022-07-21 第二版: 群友建议，1)可以用BUSCO或者RNA-seq mapping的方式来评估重复序列屏蔽的效果。 2)用水稻人工构建的高质量TE数据库评估。

最近又折腾了下植物基因组的重复序列注释，这里做一点总结


1. 我整理了许多文献（主要是植物）中关于TE占比的描述和对应的基因组大小。我发现，一般100-200Mb物种占比大概是30%以内（如拟南芥和A.lyrata)，400Mb大概是40%(如水稻）,  超过800Mb的物种占比50%以上。
2. 我对比了RepeatModeler + RepeatMasker和EDTA + RepeatMakser 这两种分析组合的结果，前者注释的比例都会比后者高。以前我之前发现一个基因组，大小也就是200Mb左右，但是TE注释比例在45.78%，。文章描述方法中给的是RepeatModeler + RepeatProteinMask/RepeatMasker的组合。我用EDTA + RepeatMakser是40.53%。
3. 然而在没有TEclass注释的前提下，RepeatModeler的TE库绝大部分注释的都是Unknown，因此RepeatMasker注释的结果中绝大部分的TE也是Unknown。我不确定这是不是因为我没有正确设置参数设置的原因。
4. 基于上一条，RepeatModeler的输出结果需要注释TE才能用于RepeatMasker。我尝试用TEclass和TEsorter对这些Unknown进行注释，两者结果相差巨大，可能1000条未知序列，TEclass会注释出900条，TEsorter只能注释出90条。
5. 我在EDTA构建的文库的基础上，分别用TEclass和TEsorter进行注释，差距会稍微小一点，但是依旧是TEclass注释的结果多。同时，对比TEclass的注释结果，和EDTA原先的注释结果，两者结果在大类上基本一致。
6. 通过查阅TEsorter的文献，发现benchmark中的sensitivity这一项的确是低于其他软件，rice的other TE只有 **16%**. 但是TEsorter的precision则优于其他软件，几乎都是100%。也就是，如果你想结果更准，TEsorter是你的不二选择，但是你想要结果更多，TEsorter就未必满足你的需求了。
7. 我发现EDTA在[F.hispida](https://ngdc.cncb.ac.cn/gwh/Assembly/7807/show)会因为TIR-finder找不到结果导致无法正确构建TE文库，不过好消息，在我自己的测序物种里面，表现都挺好的。

对于我而言，以后在重复序列注释这一块，我采取的就是EDTA + RepeatMasker这个组合。

由于我对重复序列的研究比较少，以上观点仅仅是这几天跑代码的一点心得，抛转引玉，供大家参考。