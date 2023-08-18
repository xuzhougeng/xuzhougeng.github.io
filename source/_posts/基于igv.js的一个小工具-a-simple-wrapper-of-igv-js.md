---
title: 基于igv.js的一个小工具
date: 2020-11-02 22:05:02.886
updated: 2020-11-02 22:07:39.721
url: /archives/a-simple-wrapper-of-igv-js
categories: Python
tags: 
---


当我想查看软件运行结束后得到的BAM, BigWig, BED和GTF文件时，我都需要先把他们下载到本地，然后用IGV打开。每每这个时候，我就会非常痛苦，因为我懒得下载。

为了解决这个问题，我基于igv.js 写了一个Python脚本，可以根据提供的参考基因组，bam文件输出一个网页，然后通过利用 node.js 或者 Python的网页服务器打开一个端口，直接在网页上进行查看。

代码地址: [https://github.com/xuzhougeng/myscripts/blob/master/igv_web.py](https://github.com/xuzhougeng/myscripts/blob/master/igv_web.py)

可以下面这行代码进行下载脚本, 例如我放在了家目录

```bash
wget https://raw.githubusercontent.com/xuzhougeng/myscripts/master/igv_web.py
```

```bash
python3 ~/igv_web.py  -r ref/genome.fa -m 01-read-align/*.bam -b feature.bed
```

接着用python的 SimpleHTTPServer 模块启动一个网页服务器

```bash
/usr/bin/python -m SimpleHTTPServer 9999
```

当然更加推荐使用npm的http-server，通过网页进行访问

```bash
npx http-server ./ -p 9999 -a 0.0.0.0
```

效果如下

![image.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/11/image-653ecf688d544364bc6498c0725ab215.png)
