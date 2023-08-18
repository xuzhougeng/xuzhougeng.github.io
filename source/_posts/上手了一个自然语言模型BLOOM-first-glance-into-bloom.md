---
title: 上手了一个自然语言模型BLOOM
date: 2022-06-27 10:12:07.11
updated: 2022-06-27 10:12:07.11
url: /archives/first-glance-into-bloom
categories: 数据科学
tags: 自然语言 | 深度学习
---

在2年前，OpenAI搞了一个1750亿个参数的神经网络模型, 即GPT-3(是它的前辈GPT-2, 约15亿个参数, 的100多倍)。可惜的是, GPT-3并没有开源它预训练的参数数据集，使用者只能调用他们提供的API，由于这个API必须要申请，所以我还没机会用到。不过，已经有很多公司基于GPT-3开发了一些应用，比如说GitHub就开放一个GitHub Copilot用于智能补全代码，大家感兴趣的可以去试试。

通常而言，GPT-3这种级别的神经网络模型只有OpenAI, 谷歌这种巨头科技公司能够开发和训练，毕竟1750个参数可不是一个小数目。但是，一个由约1000多个学术志愿者组成的国际志愿者不信邪，他们正在用价值700万美元的公共资金赞助的计算时间去训练一个具有1760亿个参数的自然语言模型，BLOOM。

出于兴趣爱好，我在自己的Macbook Pro上测试了这个BLOOM-1b3模型，顾名思义就是该模型有13亿参数。为什么不测试完整的1750亿个参数呢？因为光读取13亿参数到Python中，就需要占用6G左右的内存，完整的1750亿个参数，起码得是一个1T内存的服务器才行。而且完整参数版还没有完工，仍在训练中。

BLOOM是基于Transformers，而Transformers的安装要求满足如下三点

- Python 3.6+
- Flax 0.3.2+/ PyTorch 1.3.1+ / TensorFlow 2.3+ 任一框架
- Rust(用于编译tokenizers)

在Mac上，我们可以通过homebrew来安装rust

```bash
brew install rust
```

在深度学习框架上，我选择了PyTorch, 因为近期支持使用M1芯片的GPU进行加速

```bash
# 安装虚拟环境
python3 -m pip install --user --upgrade pip
python3 -m pip install --user virtualenv
# 创建虚拟环境
python3 -m venv pytorch
source pytorch/bin/activate
# 安装pytorch
pip3 install --pre torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/nightly/cpu
# 安装transformers
pip3 install transformers
```

之后打开python, 加载我们的模型(第一次运行时, transformers会自动帮忙下载模型)


```python
from transformers import AutoTokenizer, AutoModelForCausalLM

# 分词
tokenizer = AutoTokenizer.from_pretrained('bigscience/bloom-1b3')
# 模型
model = AutoModelForCausalLM.from_pretrained('bigscience/bloom-1b3')
```

由于我对自然语义处理一窍不通，所以只能根据示例代码，使用transformres的pipeline, 构建一个文本生成工具。

```python
from transformers import pipeline
generator = pipeline(task="text-generation", model=model, tokenizer=tokenizer)

# 英文输入
generator("bioinformatics is ", max_length=30, num_return_sequence=5)

# 输出结果
[{'generated_text': 'bioinformatics is  a branch of computer science that deals with the analysis of data and the design of algorithms that can be used to solve problems in'}]

# 中文输入
generator("生物信息学是一门")
# 输出结果
[{'generated_text': '生物信息学是一门新兴的交叉学科，它涉及生物信息学、计算机科学、数学'}]
```

从结果来看，13亿参数的BLOOM模型的表现已经非常符合我的预期，远超出之前摸索Transformers时测试的默认GPT-2模型。

可惜我的能力有限，对BLOOM的探索就止步于此了，等到后续我对Transformers有更多的掌握后，再更新相关的内容吧。

参考资料

- https://www.nature.com/articles/d41586-022-01705-z

