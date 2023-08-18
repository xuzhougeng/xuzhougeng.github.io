---
title: Arm架构的macOS安装conda
date: 2022-07-03 04:02:49.141
updated: 2022-07-03 04:02:49.141
url: /archives/arm64-mac-install-conda
categories: 生信软件工具箱
tags: MacOS
---

在arm架构的macOS系统安装conda实际上并没有特别需要注意的地方，我们只需要在最初下载的时候，选择conda-forge提供的miniconda，即[Miniforge](https://github.com/conda-forge/miniforge)

Miniforge是针对conda-forge优化的conda，做了如下的预设

- conda-forge作为默认channel
- PyPy的可选支持，替代标准Python(CPython)
- Mambda可选支持，替代conda
- 多种CPU架构支持(x86_64, ppc64le, aarch64(M1))


安装过程代码如下

```bash
# 下载
cd ~/Download
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
# 安装
chmod +x ~/Downloads/Miniforge3-MacOSX-arm64.sh
sh ~/Downloads/Miniforge3-MacOSX-arm64.sh
```

运行时会有一些提示信息，输入yes或者回车就行，完成后会有如下信息

![安装信息](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2022/07/image-623bd5c40e1b41e487aaddbe63c96e80.png)

注意我在途中标识的部分，这里有两个信息我们需要指导

- 我们需要重启终端才能调用conda
- conda之后会默认启动base环节，除非你用 `conda config --set auto_activate_base false` 取消这个行为。

不过，遗憾的是许多生信软件还没有arm64版本，例如bwa, samtools 等，但是没有关系这些软件可以通过 homebrew 进行安装。
