---
title: 「macOS」为什么我的时间机器的备份磁盘无法被推出？
date: 2022-06-14 09:53:24.508
updated: 2022-06-14 09:53:24.508
url: /archives/macos-why-my-driver-can-not-eject
categories: 其他
tags: MacOS
---


今天我用时间机器备份我的Mac上的数据后，打算推出我的磁盘，但是遇到了如下提示

![image.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2022/06/image-eb51acc2729b40a4bbd6183b11120c1f.png)

一般这种情况，我们有一个国际通用的解决方案，那就是关机。但是，我不想关机啊，我想找到到底是什么程序占用了这个磁盘，但是它不能顺利的被推出。

作为一个略懂linux的半个码农，我熟练的使用lsof 查看了磁盘的占用进程

```bash
$ sudo lsof  /Volumes/Time-Machine
COMMAND PID USER   FD   TYPE DEVICE SIZE/OFF NODE NAME
mds     334 root   27r   DIR   1,29      160    2 /Volumes/Time-Machine
```

结果显示，有一个PID为334的mds的进程正在使用硬盘。那么问题来了，这个mds是啥呢？ 使用百度加上 `man mds`, 我发现这是MacOS上的一个索引工具，同时我还找到一位跟我有类似经历的[文章](https://apple.stackexchange.com/questions/438329/why-is-mds-accessing-my-time-machine-drive-all-the-time)。简单说，就是因为时间机器备份需要进行文件的比对，为了实现这种快速的文件件对比，那么我们就需要建立文件索引。

强制干掉mds进程会影响我们索引的完整性，所以，最后我还是关机了。
