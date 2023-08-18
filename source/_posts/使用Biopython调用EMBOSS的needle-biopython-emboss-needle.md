---
title: 使用Biopython调用EMBOSS的needle
date: 2021-04-14 20:32:13.198
updated: 2021-04-14 20:51:03.592
url: /archives/biopython-emboss-needle
categories: 生信软件工具箱
tags: 小技巧
---

最近有一个需求，计算250条序列两两之间的相似度，之后根据相似度继续进行序列过滤，即A和B如果相似度高于一个阈值，那么就把B去掉。

我最初想到的就是在Python里面调用EMBOSS的needle，准备使用Python的os或者subprocess模块，但是突然想到之前看WGDI源代码时，发现作者时通过biopython的API进行MAFFT的调用，于是就决定改用Biopython.

通过查找文档，我找到了Biopython对这部分内容的介绍，见 <http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec97>,示例代码如下

```python
from Bio.Emboss.Applications import NeedleCommandline
needle_cline = NeedleCommandline(asequence="alpha.faa", 
                                 bsequence="beta.faa", 
                                 gapopen=10, gapextend=0.5, 
                                 outfile="needle.txt")
stdout, stderr = needle_cline()
align = AlignIO.read("needle.txt", "emboss")
```

但是这个代码距离我最终代码还有很长的路要走，因为我需要先将序列写出到文件，然后调用命令行，然后再读取输出文件。我的问题是，能不能省掉这几次I/O操作。在Emboss文档的下面，我发现Biopython提供了两个模块，分别是Bio.pairwise2和Bio.Align.PairwiseAligner, 可以直接使用Biopython的.SeqRecord对象。按照文档描述，pairwise2和PairwiseAligner能够达到和调用EMBOSS相同的速度，经过我测试，发现两者速度的确差别不大，这说明直接在Python里面进行配对序列联配并不比先输出运行所需文件，然后读取运行生成文件这个流程快，所以经过一番折腾之后，我还是决定使用NeedleCommandline。

为了方便后续调用，我写了一个函数，函数的前3个参数是用来从SeqRecord列表中提取序列用于写出，第4个参数用来提供needle程序的路径。

```python
def run_needle(aseq_pos, bseq_pos, seq, needle_path, verbose=False):
    """run EMBOSS needle retrun align 
    """

    # get the Bio object
    aseq = seq[aseq_pos]
    bseq = seq[bseq_pos]

    # create tempfile for write fasta and needle result
    aseq_file = mkstemp(suffix=".fas")[1]
    bseq_file =  mkstemp(suffix=".fas")[1]

    needle_file =  mkstemp(suffix=".needle")[1]

    # write the sequence
    SeqIO.write(aseq, aseq_file, "fasta")
    SeqIO.write(bseq, bseq_file, "fasta")

    # build needle commnad
    if needle_path is None:
        needle_cmd = NeedleCommandline( asequence=aseq_file, bsequence=bseq_file,
                   gapopen=10, gapextend=0.5, outfile=needle_file)
    else:
        assert os.path.isfile(needle_path)
        needle_cmd = NeedleCommandline( needle_path, asequence=aseq_file, bsequence=bseq_file,
                   gapopen=10, gapextend=0.5, outfile=needle_file)
    
    # running needle_cmd
    stdout, stderr = needle_cmd()
    if verbose:
        print(stdout + stderr, file = sys.stderr)
    
    align = AlignIO.read(needle_file, "emboss")
    # deleting the tempfile
    os.unlink(aseq_file)
    os.unlink(bseq_file)
    os.unlink(needle_file)
    return align
```

在写这个函数时，我在输出文件的文件名上折腾了很长时间，最初想使用原序列的ID，但ID包括一些奇奇怪怪的字符，例如`(|:`, 而NeedleCommandline不会自动在构建的命令中增加引号，使得运行时会出错。(你可以试试`echo a(b`会出现什么情况 )

后来，我使用了随机数作为文件名。这个方法在单个进程时问题不大，当时当我为了加速使用多进程时，却发现发现随机数重复了。进程A输出的文件A由于和进程B的输入文件名一样，导致被跑完的进程B删了。当然，解决方法也很简单，加上一些判断文件是否存在的ifelse语句即可。

不过，我还是选择调用tempfile标准库，使用mkstemp获取临时文件名，才解决了输出文件名的危机。

最终我就可以在循环中调用这个函数，然后将结果保存在列表中，用于后续处理。

```python
fasta = "/path/to/your/sequence.fasta"
seqs = [ fas for fas in SeqIO.parse(fasta, "fasta") ]
needle_path = "/opt/biosoft/EMBOSS-6.6.0/bin/needle"

align_res = []

for i in range(len(seqs)):
    for j in range(i+1, len(seqs)):
        align = run_needle(i, j, seqs , needle_path)
        align_res.append(align )
# 以pickle保存python对象
outfile = open("pairwise_stats.txt", "w")
for align in align_res:
    line = get_psa_stat(align)
    line = "\t".join(line) + "\n"
    outfile.write(line)

outfile.close()
```

根据我的计算，循环要运行3万多次，按照平均1s一个循环，只要500多分钟，小于9个小时，我只要在回去前运行这个程序（写完第一版大概是晚上10点），然后第二天早上过来收结果就行了。

但是，正如前文提及到的，「该程序会在多个进程运行时报错」，你也就知道我没有跑循环，而是选择了通过并行对循环进行加速。在一篇中，我将会讲一个，如何把一个简单的问题复杂化，最终又回到最初代码的故事


