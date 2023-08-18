---
title: 如何高效地从BAM文件中提取fastq
date: 2019-07-29 10:05:32.347
updated: 2019-08-24 18:03:41.058
url: /archives/如何高效地从BAM文件中提取fastq
categories: 生信软件工具箱
tags: 
---

在[如何从BAM文件中提取fastq](https://www.jianshu.com/p/368b9471656e)文章里其实也发现了从BAM里面提取Fastq是有些麻烦，只不过最后通过`samtools`的子命令实现了数据提取，实现功能之后也没有再去思考如何提高效率。

最近读到[每周文献-190419-植物单细胞BAM重比对以及假基因研究](https://www.jianshu.com/p/e866ae780e79)时，发现里面提到了一个工具叫做 **bazam**, 功能就是提取Fastq文件，文章发表在 Genome Biology 上。

软件地址为[https://github.com/ssadedin/bazam](https://github.com/ssadedin/bazam)，因为他是个Groovy工具（据说Groovy比Java好写）所以安装很方便，

```bash
git clone https://github.com/ssadedin/bazam.git
cd bazam
git submodule update --init --recursive
./gradlew clean jar
```

可以通过如下命令检查是否安装成功

```bash
java -jar build/libs/bazam.jar -bam  test.bam  > tmp.fastq
```

假如你原来的BAM文件是双端测序，想提取成两个文件，那么可以用如下命令。**注**：这里的your.bam 指的是你实际BAM文件

```bash
java -Xmx16G -jar build/libs/bazam.jar \
    -bam your.bam -r1 your_R1.fq -r2 your_R2.fq &
```

此外，bazam还支持使用`-L`指定区间来提取某个区间的数据。

> 尽管是可以从BAM文件中提取Fastq文件，但如果原始BAM文件过滤了未必对的read，那么你还是会损失掉这部分信息，建议保存好原来的Fastq文件。

这里补充下为啥从BAM文件中提取双端序列那么麻烦，这是因为一般而言BAM文件都是按照位置信息排序，想要找到配对的reads，要么是根据read的编号进行排序（这个方法要求额外的内存和存储空间），或者就是在提取的时候记录当前的read的ID，再找到另一端ID后释放内存空间。

## 参考资料

- [Bazam: a rapid method for read extraction and realignment of high-throughput sequencing data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1688-1)
- [https://github.com/ssadedin/bazam](https://github.com/ssadedin/bazam)
