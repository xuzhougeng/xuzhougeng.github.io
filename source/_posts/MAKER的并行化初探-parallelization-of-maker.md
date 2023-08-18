---
title: MAKER的并行化初探
date: 2020-06-29 04:28:13.909
updated: 2020-07-07 07:16:18.851
url: /archives/parallelization-of-maker
categories: 基因组学
tags: 注释 | 流程工具 | MAKER
---


MAKER并行分为两种，一种基于MPI，运行方式为`mpiexec -n 线程数 maker`, 一种是在同一个项目中运行多次maker。前者需要在安装MAKER时进行设置，后者相当于你手动按照染色体数目进行拆分，然后分开运行MAKER。本片主要介绍基于MPI的并行策略。

下图是MAKER是基于MPI的并行化流程。分为三个层次，contig水平，更小的片段水平，不同程序的并行。

contig水平就是对每条contig进行分别注释。而由于每个contig的长度不一，当较短的contig运行结束后，我们还需要等待较长的contig结束。因此，我们还需要将contig继续拆分，让每条contig分成更小的片段，成为基本的分析单元。对于每个分析单元，还可以考虑将不同程序进行并行。

![Fig1](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/06/image-5624eadd6cc346a4b93a63be53fea0fb.png)

在MAKER实际运行时，如果用gotop进行查看你会发现系统中会出现很多maker进程(以`mpiexec -n 10`为例)

![Fig2](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/06/image-004af64038c349fb99dbb3ee58531b58.png)

但是你发现只有少量的和注释有关的程序(如blastx)。如果用`ps aux | grep maker`查找和maker有关的进程，你会发现大量的`maintain.pl`进程，并且这些进程后面跟着一大串你不认识的符号

```bash
/usr/bin/perl /opt/biosoft/maker/bin/../lib/File/maintain.pl 46535 30 %04%09%08%31%3...
```

为了理解MAKER的并行，你需要通过`pstree -ap [mpiexec进程ID]`更细致的了解MAKER的进程树

![Fig3](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/06/image-e6e682747a3a443e9a2b0e9810f0177c.png)

从中你可以发现，因为我们用`mpiexec`启动了10个maker进程，所以`hydra_pmi_proxy`有10个maker子进程。对于这10个maker子进程，每个进程都至少会有一个`maintain.pl`。通过阅读`maintain.pl`的源代码，我们可以得知该程序后接参数中的46535是PID(进程号), 30是sleep time, 最后一个是URI编码字符，是利用`Storable::freeze`持久化的数据结构，可以通过如下的代码进行解码（命名为decode.pl)

```perl
#!/usr/bin/env perl
use warnings;
use strict;

use URI::Escape;
use Storable;
use vars qw($LOCK);

my $serial = shift;

$serial = uri_unescape($serial);
$LOCK = Storable::thaw($serial);

print "$LOCK->{lock_file} \n"
```

可以批量对`maintain.pl`里的信息进行解码

```bash
ps aux | grep maintain.pl | grep -v grep | awk '{print $15}' | xargs -i perl decode.pl {}
```

因此，那么多的maker并不是在话说，而是在负责任务调度，差不多一个任务会有3个maker进程保驾护航。

回到之前的进程数，这里有一些maker进程下会跟着一些正在运行的命令，blastx, exonerate

![Fig4](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/06/image-01a3786169a74d278d804a79c9e36e91.png)

而有些只有maker和maintain.pl进程。

![Fig5](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/06/image-6d38421ce79d4c129fe92e78d71a926d.png)

如果等待一会，再次运行`pstree -ap [mpiexec进程ID]`，你会发现这两种状态会发生转换，这说明任务启动也需要一段时间。
