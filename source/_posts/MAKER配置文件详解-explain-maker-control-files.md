---
title: MAKER配置文件详解
date: 2020-07-01 06:26:12.241
updated: 2020-07-07 07:16:03.292
url: /archives/explain-maker-control-files
categories: 基因组学
tags: 注释 | MAKER
---


本文翻译自<http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/The_MAKER_control_files_explained>

MAKER会生成三个配置文件，如下

- maker_opts.ctl: 控制MAKER分析行为的主要配置文件
- maker_bopts.ctl:  BLAST和Exnerate的过滤阈值
- maker_exe.ctl: MAKER运行过程中所依赖软件的路径

## maker_opts.ctl

### 基因组

用于设置被注释的基因组序列的位置和物种类型，包括`genome`和`organisam_type`两项

- genome: FASTA序列路径
- organism_type: eukaryotic,prokaryotic二选一

需要注意的是，基因组序列的N50需要超过预期基因长度的中位数，否则注释效果不好。另外最好保证基因组序列只包括A,T,C,G,N, 对于其他类型兼并碱基可以都改成N.

### 使用MAKER得到GFF3进行重注释

这一项基本上我们用不上，它是在当你把MAKER的中间输出文件都删除了，仅保留了输出的GFF3文件时，你可以用之前相同的输入设置重新运行流程得到相同的输出。

```bash
#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no
```

### EST/转录本证据

出于历史原因, MAKER还是用EST代表了之前的EST数据和目前的转录组数据。 此处不只是使用EST数据，而是可以使用组装的mRNA-seq, 组装的全长cDNA。 我们预期他们能够正确的组装，并联配到正确的剪切位点(对于FASTA格式，MAKER使用`exonerate`找到剪切位点）。用途如下:

- 直接**推断**基因模型(est2genome)
- 作为预测结果的支持证据
- 修改结果和增加UTR
- 鉴定可变转录本
- 在某些情况下，这些数据和其他证据能帮助MAKER推断基因边界
- 在预测步骤中辅助基因预测工具推断剪切位点

设置项如下:

- `est`: EST, 全长cDNA, 组装的mRNA序列，可以通过逗号分隔多个文件。
- `altest`: 如果真的没有当前物种的转录组数据，也可以使用同源物种的序列。这些序列会通过tblstx进行比对，会消耗大量的运算时间。
- `est_gff`: 预比对的转录本，以GFF3格式存放，通常来自于CuffLinks/StringTie的组装结果，或者是上一步的MAKER注释
- `altest_gff`: 和altest一样，只不过是比对后以GFF3存放，基本上也用不到

### 同源蛋白证据

和之前的转录组数据类似，用途如下

- 直接推断基因模型(protein2genome), 仅在它们能够正确的联配到剪切位点附件。
- 作为预测结果的支持证据(MAKER会检查CDS，保证基因预测结果和蛋白联配是相同的阅读框)
- 某些情况下，用于推断基因边界
- 在预测步骤中，使用从蛋白推断的ORF辅助从头预测软件

**建议**使用 uniprot/swiss-prot 或 RefSeq上的NP数据，因为经过人工审查，可信度较高。不建议是用UniProt/tremble或者Genbank上的数据，这些数据的可信度较低。你可以挑选几个同源物种的高可信度蛋白。或者使用MAKER注释的其他物种AED小于0.5的转录本产物。

由于许多注释里包含一些死亡转座子(dead transposons)或伪基因(pseudogenes)，因此**不建议**使用临近物种的所有注释蛋白。我们想象一个比较糟糕的情况，如果你有邻近物种的死亡转座子，当你构建你的重复序列屏蔽文库时，你发现其中一个条目和该序列匹配。 于是你假设这是一个真实的基因，于是你从屏蔽文库中删除了该条目。吸纳子啊，当你注释基因组的时候，该基因变成了注释集中的一整个基因组家族，但这其实是糟糕的证据和重复序列屏蔽所导致的后果。

需要设置的就两项：

- `protein`: 以FASTA存放的蛋白序列位置
- `protein_gff`: 以GFF3格式记录的蛋白序列预先比对结果，通常来自于之前MAKER输出

### 重复序列重复屏蔽

我们可以通过屏蔽重复序列来避免EST和蛋白比对到重复区域，防止基因预测算法在这些区域预测外显子。由于许多重复序列会编码真实的蛋白(例如反转座子等)，基因预测工具和比对工具会被他们所迷惑(会在一个基因中错误的加上外显子)

- `model_org`: 默认是all, 可以设置成RepeatMasker中RepBase数据库的其中一个
- `rmlib`: 自己构建的重复序列文库
- `repeat_protein`: 转座因子的蛋白序列
- `rm_gff`: 以GFF3格式记录的重复序列位置信息
- `prok_rm`: 不需要修改，因为原核生物不需要考虑重复序列
- `softmask`: 不需要修改

### 从头基因预测

如果你需要从MAKER以外获取基因模型，则需要在这一节添加相应的配置。根据可信度高低，MAKER会对这些基因模型采取不同的行为。

- 通过软件预测的基因结构可信度低，它们不会影响证据簇(evidence cluster). MAKER会保留预测结果，或者根据EST证据调整外显子，如果有证据支持，那么他们会保留在最终的注释集中。如果有多个注释结果，MAKER会对其进行比较，从中挑选出最优结果。
- `model_gff`提供的基因模型的可信度最高，会影响证据簇。在一些基因边界判定中，MAKER在证据簇的影响下，更倾向于保留之前的基因模型而非替换。它们也会保留名字。MAKER不会修改模型，要么删除要么保留。最后，即便没有证据支持，MAKER还是会保留他们，而不会删除他们，只不过最终的AED会设置为1.

- snaphmm: SNAP的HMM文件路径，允许多个输入，以逗号分隔
- gmhmm: GeneMark的HMM文件文件路径，允许多个输入，以逗号分隔
- augustus_species: Augustus的基因模型命名, 不容易训练，但是效果很好
- fgenesh_par_file: FGENESH的HMM参数文件，收费工具，基本上用不到
- pred_gff: 其他预测工具的输出结果，以GFF3格式保存
- model_gff: 最高可信度的GFF输入
- run_evm=0: 是否让MAKER运行EVM，速度会变慢
- est2genome: 让MAKER根据EST推测基因模型
- protein2genome: 让MAKER根据蛋白序列推测基因模型
- trna: 使用tRNAscan分析tRNA #find tRNAs with tRNAscan, 1 = yes, 0 = no
- snoscan_rrna: Snoscan分析snoRNA所需的rRNA文件= #rRNA file to have Snoscan find snoRNAs
- snoscan_meth：Snoscan分析snoRNA所需的O-methylation site 文件
- unmask: 是否在屏蔽序列中进行基因预测，默认是0
- allow_overlap: 允许的基因重叠比例，从0到1，空白表示使用默认值

### 其他类型的注释

这一项功能很简单，就是提供一个GFF文件，在MAKER运行结束后增加里面的信息

- `other_gff`: 其他类型注释的GFF文件路径

### 外部程序选项

这里的两个参数用于影响外部程序，即BLAST的行为

- `alt_peptide`: 对于非标准氨基酸的替换方法，默认是C(cysteine)
- `cpus`: BLAST的线程数，如果使用MPI，该值可以设置的小一些，默认是1.

### MAKER行为选项

这里的选项用于调整MAKER的行为，使其符合你的基因组特性

- `mad_dna_len`至少要3倍于预期的最大内含子长度。在内存足够的情况下，对于**脊椎动物**可以考虑设置为`300000`，植物一般没有那么大的内含子
- `min_contig`: 低于该值的contig会被过滤掉，建议设置为10k.
- `pred_flank`: 在基因预测时，将证据簇在两端进行扩展，默认是200 bp. 对于比较紧凑的基因组，降低该值能够避免基因错误合并。对于比较稀疏的基因组，提高该值可以避免外显子缺失。
- `pred_stats`: 默认是0，只计算MAKER预测基因的AED和QI值，设置为1则计算所有从头预测的基因结构。
- `AED_threashold`:  根据AED(0-1)值来过滤输出的基因，默认是1，表示保留所有预测结果。
- `min_protein`:  一些时候，基因预测工具会生成许多短预测结果，而由于一些证据类型（例如mRNA-seq）存在噪音，导致这些预测结果看起来有证据支持，于是保留在最终的输出结果中(AED>1). 通过限制预测蛋白的氨基酸数(amino acides, AA)，可以减少预测结果。
- `alt_splice`: 是否计算可变转录本
- `always_complete`: 这个是MAKER开发者在合作者的要求下加上的参数，用来确保基因模型始终有起始密码子和终止密码子，默认是0
- `map_forward`: 用于保留老版本的GFF文件的信息，映射到新的版本GFF中。
- `keep_preds`:  设置为1时表示保留所有的预测结果，默认是0.
- `split_hit`:   新版本的MAKER(>2.28)不需要考虑该项。因为之前版本拆分后的contig是互不重叠，于是就有可能有外显子被刚好被拆成两端，设置该项可以保留该信息。
- `single_exon`: 如果一个EST只有只有单个外显子，默认情况下MAKER并不会把它当做支持基因模型的证据, 除非还有同源蛋白作为支持。单外显子EST和组装的mRNA-seq转录本通常是RNA制备过程中的基因组序列污染。关闭时会降低MAKER的敏感度(sensitivity), 但是当你打开它的时候，命中的特异性会比整体的准确度(accuracy)差得多。
- `single_length`: 在启用`single_exon`时，设置该项用于保留一些比较小的序列，但即便注释了也可能不是有功能的蛋白。
- `correct_est_fusion`用来避免因为UTR的重叠导致将基因模型的错误合并，在**真菌基因组**中比较常见。它会检查基因模型的5' UTR长度是否超过基因长度的一半，如果是的话，那么MAKER会在起始密码位置打断基因，然后在5'UTR区重新预测基因。
- `tries`: 尝试次数，默认是2
- `clean_try`: 重新尝试时，是否删除之前的文件，默认是0，也就是不删除。建议设置为1
- `clean_up`: 删除分析过程中的文件，默认不删除
- `TMP`: 临时文件的目录，默认存放在`/tmp`下，建议设置一个容量比较大的目录