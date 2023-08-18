---
title: 深入探索biopython的配对序列联配模块
date: 2022-01-02 08:34:48.384
updated: 2022-01-02 08:34:48.384
url: /archives/explore-biopythonpairwise-alignment
categories: Python
tags: 算法
---


> 本文基于biopython>=1.79

配对序列联配是一个非常基础的生信分析，biopython提供了Bio.Align.PairwiseAligner基于动态规划处理两个序列的联配问题。

如下是一个简单的案例：

```python
from Bio import Align
aligner = Align.PairwiseAligner()
alignments = aligner.align("ACCGGTTT", "ACGGGTT")
```

它首先初始化了一个PairwiseAligner对象，然后用该对象对两条序列进行联配，联配结果记录在了alignments。

alignments里面记录得分最高的联配结果（不止一个），我们可以通过循环来遍历其联配结果。

```python
for alignment in alignments:
    print("Score = %.1f:" % alignment.score)
    print(alignment)
```

PairwiseAligner的运行速度极快，这是因为他底层通过C语言实现，而非Python代码，如下是他的结构体。

```c
typedef struct {
    PyObject_HEAD
    Mode mode;
    Algorithm algorithm;
    double match;
    double mismatch;
    double epsilon;
    double target_internal_open_gap_score;
    double target_internal_extend_gap_score;
    double target_left_open_gap_score;
    double target_left_extend_gap_score;
    double target_right_open_gap_score;
    double target_right_extend_gap_score;
    double query_internal_open_gap_score;
    double query_internal_extend_gap_score;
    double query_left_open_gap_score;
    double query_left_extend_gap_score;
    double query_right_open_gap_score;
    double query_right_extend_gap_score;
    PyObject* target_gap_function;
    PyObject* query_gap_function;
    Py_buffer substitution_matrix;
    PyObject* alphabet;
    int* mapping;
    int wildcard;
} Aligner;
```

在结构体中，有一个非常关键的algorithm, 决定了联配时所需的算法。主要分为三个大类 **NeedlemanWunschSmithWaterman**, **WatermanSmithBeyer**和**Gotoh**.算法无法手动指定，而是根据match, mismatch, gap等参数的得分来确定(共6种)，具体的决策过程如下

当gap得分与gap长度有关，也就是 gap_score 赋值一个**函数**时，在WatermanSmithBeyer基础上，当mode是global, 使用 "Waterman-Smith-Beyer global alignment algorithm"；当mode是local，使用 "Waterman-Smith-Beyer local alignment algorithm"

当gap打开的得分和gap延伸的得分相同时，在NeedlemanWunschSmithWaterman基础上，如果mode是global， 使用Needleman-Wunsch；如果mode是local， 使用Smith-Waterman

其他情况下，也就是得分不同，则是在Gotoh基础上，如果mode是global，选择 "Gotoh global alignment algorithm"，如果mode是local，则选择"Gotoh local alignment algorithm"


我们举几个简单的例子来说明

场景1: 通过函数定义为gap打分

```python
from Bio import Align
aligner = Align.PairwiseAligner()
aligner.gap_score = lambda x: x*2
aligner.algorithm
# 'Waterman-Smith-Beyer global alignment algorithm'
aligner.mode = "local"
aligner.algorithm
# 'Waterman-Smith-Beyer local alignment algorithm'
```

场景2: gap打开和延伸得分相同

```python
from Bio import Align
aligner = Align.PairwiseAligner()
aligner.open_gap_score = -1
aligner.extend_gap_score = -1
aligner.algorithm
# 'Needleman-Wunsch'
```

场景3: gap打开和延伸得分不同

```python
from Bio import Align
aligner = Align.PairwiseAligner()
aligner.open_gap_score = -1
aligner.extend_gap_score = -2
aligner.algorithm
# 'Gotoh global alignment algorithm'
```


除了为match和mismatch设置固定值外，我们还可以考虑使用得分矩阵为不同碱基或者氨基酸设置不同的得分，例如BLOSUM62

```python
from Bio import Align
aligner = Align.PairwiseAligner()
# add matrix 
from Bio.Align import substitution_matrices
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
# align
alignments = aligner.align("ACCGGTTT", "ACGGGTT")
```

biopython的替换矩阵来自于<ftp://ftp.ncbi.nih.gov/blast/matrices/>, 保存在Bio/Align/substitution_matrices/data 目录下


在搞定联配之后，我们如何获取联配后的结果？实际上对于一个 PairwiseAlignment 对象，例如 alignment, 我们调用 `print(alignment)`实际等价于 `print(format(alignment))`, 因此我们可以通过如下代码，获取联配后的query和target序列

```python
seqs = format(alignment).split("\n")
target = seqs[0]
query = seqs[2]
```

最后总结下本文内容：

1. biopython的配对序列联配算法是C编写，运行速度极快
2. 联配算法会根据 gap score和mode 自动确定
3. 基于format函数可以获取联配后的序列
