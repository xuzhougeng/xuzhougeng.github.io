---
title: 「小技巧」如何让IGV更快的加载GTF和GFF注释文件
date: 2021-01-23 19:03:03.85
updated: 2021-01-23 19:03:03.85
url: /archives/trick-load-gtf-and-gff-to-igv-more-quickily
categories: Linux
tags: 小技巧
---


很简单，就下面3行命令

```bash
gff=
(grep ^"#" $gff; grep -v ^"#" $gff | sort -k1,1 -k4,4n) | bgzip > sorted.gff.gz;
tabix -p gff sorted.gff.gz;
```

第一行的gff是定义输入文件。第二行是对GFF文件进行排序。第三行是利用HTSLIB中的tabix工具建立索引，得到一个sorted.gff.gz.tbi 索引文件。

后续IGV在加载文件时，会根据索引文件直接读取给定区域的信息，不但可以降低内存的占用，还可以提高加载速度。