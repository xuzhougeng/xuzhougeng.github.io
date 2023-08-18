---
title: 从一则12年前的提问中学习:从配对序列联配到多序列联配
date: 2022-10-21 13:48:24.645
updated: 2022-10-21 13:48:24.645
url: /archives/qa-from-pairwise-alignments-to-multiple-alignment
categories: 基因组学
tags: 比较基因组学
---

最近学习多基因组比对时，看到一则12年前在Biostars发布的提问, [Programming Challange: Pairwise Alignments To Multiple Alignment](https://www.biostars.org/p/2882/), 收获颇多，这里记录下。

提问者，一开始阐释了自己的问题，也就是他有10-12个非常近的物种的染色体序列，将这些物种和一个参考染色体比对后，得到了多个结果。他希望，在生成多序列联配结果的同时不影响到原本单独的比对结果。那么，他想的就是，在各个序列中加上一些插入，就可以得到相对基因组的全局联配。

同时，他强调了，自己不是来找序列相似度！他想的是有没有已有的脚本可以做这些事情，或者提供一些代码上的建议帮助他完成。

甚至，他还给了一个案例，来说明自己的需求

也就是把下面这段

```text
Ref1: CGACAAT--GCACGACAGAGGAAGCAGAACAGATATTTAGATTGCCTCTCATTTTCTCTCCC
Seq1: CGACAATAAGCACGACAGAGGAAGCAGAACAGATA-----ATTGCCTCTCATTTTC-CTCCC

Ref1: CGACAATGCACGACAGAGGAAGC--AGAACAGATATTTAGATTGCCTCTCATTTTCTCTCCC
Seq2: CGACAAT-CACGACAGAGGAAGCTTAGAACAGATATTTAG---GCCTCTCATTTTCTCTCCC

Ref1: CGACAATGCACGACAGAGGAAG----CAGAACAGATATTTAGATTGCCTCTCA----TTTTCTCTCCC
Seq3: CGACAATGCACGACAGAGGAAGTTTTCAGAACAGATATTTAGATTGCCTCTCAAAAATTTTCTCTCCC
```

变成下面这段。

```text
Ref1: CGACAAT--GCACGACAGAGGAAG----C--AGAACAGATATTTAGATTGCCTCTCA----TTTTCTCTCCC
Seq1: CGACAATAAGCACGACAGAGGAAG----C--AGAACAGATA-----ATTGCCTCTCA----TTTTC-CTCCC
Seq2: CGACAAT---CACGACAGAGGAAG----CTTAGAACAGATATTTAG---GCCTCTCA----TTTTCTCTCCC
Seq3: CGACAAT--GCACGACAGAGGAAGTTTTC--AGAACAGATATTTAGATTGCCTCTCAAAAATTTTCTCTCCC
```

我觉得大部分人看到这样子的提问，就都知道提问者到底需要什么，也就不需要花太多时间思考题目，问提问者更多细节。


排名第一的回答来自于唐海宝老师

他首先回答了，作者需要的工具是TBA/MULTIZ, 可以从[Miller Lab](http://www.bx.psu.edu/miller_lab/)下载.

接着补充了一点细节，多序列联配的潜在原则是，gap和insertion的引入和你比对序列的顺序有关。也就是说，在很多情况下，`seq1-seq2-seq3和`seq2-seq1-seq3`结果是不一样的，“once a gap, always gap”

最后说了TBA软件的不足，即需要用他们定义的MAF格式作为输入，也就是用户得做一些格式转化你工作。并给了一个使用案例

已知参考序列是ref1, 用于比对的序列是seq1, seq2, seq3。比对之后得到ref1.seq1.sing.maf, ref1.seq2.sing.maf, ref1.seq3.sing.maf, 这三个文件。提供一个进化树描述序列的顺序，如(((ref1 seq1) seq2) seq3， 表示ref1和seq1近，后面跟seq2近，最后是seq3。

最后运行如下命令

```bash
tba "(((ref1 seq1) seq2) seq3)" *.*.maf tba.maf
```

输出的tba.maf 就是你想要的结果，Good luck!


