---
title: 放弃conda拥抱mamba
date: 2022-08-03 02:22:18.863
updated: 2022-08-03 05:06:53.513
url: /archives/use-mamba-instead-of-conda
categories: Python
tags: conda
---


最早conda是Anaconda提供的包管理工具，用于管理python相关库，以及安装一些和Python相关的C/C++环境。

既然C/C++环境能够配置，那么Java, Fortran, JavaScript等编程语言应该也可以支持，于是conda能够管理的软件包变得越发的多。尤其是社区引导的开源项目[conda-forge](https://conda-forge.org)的出现，几乎市面上所有开源工具都可以通过conda进行管理。

但一个可怕的问题也随之而来，那就是原本的包管理工具conda执行速度实在是太慢了，无法处理那么庞大的依赖环境。尤其是当我想去安装一个生物信息学的分析流程时(`conda create -c conda-forge -c bioconda -n EDTA EDTA`)，conda的安装界面转了一个晚上都没能处理完。为了处理这一危机，mamba孕育而生。它始于2019年，由Wolf Vollprecht开发，拥有和conda完全一样的使用体验，你完全可以用 `alias conda=mamba` 让mamba作为conda的别名。虽然使用方法完全一样，但是在速度上，mamba则是远超自己的前辈，还是以之前的EDTA安装为例，mamba只需要不到1分钟就可以完成依赖环境的解析，同时它还支持异步下载安装，使得整体速度都有一个质的飞跃。

那么问题来了，如何才能安装配置mamba呢？这里有两个方法。

当你已经安装了miniconda或者anaconda时，我们可以通过conda管理工具去安装, 即 `conda install -c conda-forge -n base mamba`

如果是在一台全新的服务器上准备安装conda，那么我建议直接配置miniforge作为miniconda的替代。例如下面的代码用于安装mamba作为包管理工具, PyPy3.7作为解释器conda环境。

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-pypy3-Linux-x86_64.sh
bash Mambaforge-pypy3-Linux-x86_64.sh
```


为什么要用miniforge而不用原来的miniconda呢？主要原因是miniforge会默认配置好conda-forge这个channel, 有些生物信息学软件依赖的包来自于conda-forge而非默认channel。如果没有这个默认配置，可能安装的软件还不能用，例如samtools.