---
title: 吐槽下biopython的SeqIO.write导出fasta
date: 2022-02-17 08:40:38.045
updated: 2022-02-17 08:40:38.045
url: /archives/biopython-seqio-write-are-not-good
categories: Python
tags: 流程工具
---

当使用biopython的SeqIO模块加载fasta之后，每条fasta都会记录成 `Bio.SeqRecord.SeqRecord`, 会有id, name, description, dbxrefs 等属性.

举个例子，输入数据是 test.fa

```python
# test.fa
#>a|b c
#ATCG
seq = SeqIO.read("test.fa", "fasta")
seq
SeqRecord(seq=Seq('ATCG'), id='a|b', name='a|b', description='a|b c', dbxrefs=[])
```

我们可以将其保存成另一个文件

```python
with open("test2.fa", "w") as f:
    SeqIO.write(seq, f, "fasta")
```


输出的header和输入的header是完全一样的。


假如，你想修改他的名字, 比如说把 "a|b" 改成 "a", 你会发现一个非常诡异的事情，就是它完全不会按照你想的来.

```python
seq.id = "a"
# 输出 a a|b c
seq.name = "a"
# 输出a|b c
seq.description = "a c"
# 输出 a|b a c
```

如果你改了他的ID，他会输出 ID + description, 如果你改了name，相当于没改，如果你改了description，就真的只有description（即header中空格后的信息）会发生变化。

假如我想输出 ">a c", 那我就必须下面这个形式

```python
seq.id = "a"
seq.description = "c"
```

我实在是没想懂他背后的逻辑, 修改输出的header居然如此不符合预期
