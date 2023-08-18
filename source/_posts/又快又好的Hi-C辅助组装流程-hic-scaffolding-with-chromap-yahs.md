---
title: 又快又好的Hi-C辅助组装流程
date: 2023-01-27 02:38:28.742
updated: 2023-01-27 02:50:31.542
url: /archives/hic-scaffolding-with-chromap-yahs
categories: 生信软件工具箱
tags: 基因组
---


之前使用的是3D-DNA流程做Hi-C的辅助组装，它的最大优势就是输出结果可以对接下游的JBAT(juicerbox with Assembly Tools)进行手动矫正。然而它点缺陷也很明显，处理速度不够快，且对植物的优化不行，同时目前许久不更新了。

最近我发现了一套新的组合，chromap + yahs 完全替代之前3D-DNA流程。它的依赖工具如下

- chromap: 高效Hi-C数据比对
- samtools: sam转bam
- yahs: 另一个Hi-C scaffolding工具。纠错上准确性高，排序上略强3d-dna，远超SALSA2。
- juicer_tools: 用于输出导入JuiceBox

chrompa, samtools, yahs可以直接用conda进行安装，juicer_tools依赖Java环境，并需要单独下载

```
conda create -n hic-scaffolding -c bioconda -c conda-forge chromap samtools yahs assembly-stats openjdk 
# 1.19.02版本就行了, 最新的3.0不向下兼容
wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.19.02.jar
```


具体分析步骤如下，我们需要提供前期组装结果，以及Hi-C数据

```bash
contigsFasta=/到/你的/contig.fa的路径
r1Reads=/到/你的/Hi-C R1测序的路径
r2Reads=/到/你的/Hi-C R2测序的路径
```


第一步，数据比对

```Bash
# 建立索引
samtools faidx $contigsFasta
chromap -i -r $contigsFasta -o contigs.index

# alignment
chromap \
    --preset hic \
    -r $contigsFasta \
    -x contigs.index \
    --remove-pcr-duplicates \
    -1 $r1Reads \
    -2 $r2Reads \
    --SAM \
    -o aligned.sam \
    -t 50

# 排序   
samtools view -bh aligned.sam | samtools sort -@ 50 -n > aligned.bam
rm aligned.sam    
```

> 按照read的名字进行排序和按照位置排序或未排序的结果会有一些不同。

第二步，又快又好的scaffolding。默认只需要两个输入，组转的contig.fa和比对的bam，和C语言一样简洁。

```Bash
yahs $contigsFasta aligned.bam
```

在输出结果中

- inital_break 表示纠错的中间输出
- _scaffolds_final.agp和_scaffolds_final.fa则是最终结果

对于输出结果，我们希望进行可视化，此时可以使用yahs提供的jucier工具

第三步，为juicer_tools准备输入

```Bash
juicer pre -a -o out_JBAT \
    yahs.out.bin \
    yahs.out_scaffolds_final.agp \
    $contigsFasta.fai
# -o out_JBAT表示输出文件名的前缀    
```

一共包括如下几个文件

- out_JBAT.assembly
- out_JBAT.assembly.agp
- out_JBAT.hic
- out_JBAT.liftover.agp
- out_JBAT.txt

out_JBAT.txt就作为下游的输入

```Bash
JUICER=/路径/到/juicer_tools_1.19.02.jar
asm_size=$(awk '{s+=$2} END{print s}' $contigsFasta.fai)
java -Xmx36G -jar $JUICER \
    pre out_JBAT.txt out_JBAT.hic <(echo "assembly ${asm_size}")
```

输出的out_JBAT.hic就可以导入到JBAT进行组装,导出为out_JBAT.review.assembly 

将手动修改的结果传递给juicer，进行scaffolding。

```Bash
juicer post -o out_JBAT out_JBAT.review.assembly out_JBAT.liftover.agp contigs.fa
```

输出结果为 out_JBAT.FINAL.agp, out_JBAT.FINAL.fa
