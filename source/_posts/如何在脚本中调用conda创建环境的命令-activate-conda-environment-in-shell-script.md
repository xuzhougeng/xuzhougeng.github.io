---
title: 如何在脚本中调用conda创建环境的命令
date: 2022-11-21 01:04:58.225
updated: 2022-11-21 01:04:58.225
url: /archives/activate-conda-environment-in-shell-script
categories: Linux
tags: conda
---

事情源于，我在写脚本的时候，在脚本里面插入了一句 `conda activate 环境名`, 然后出现如下的报错提示

```text
CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'.
To initialize your shell, run

    $ conda init <SHELL_NAME>

Currently supported shells are:
  - bash
  - fish
  - tcsh
  - xonsh
  - zsh
  - powershell

See 'conda init --help' for more information and options.

IMPORTANT: You may need to close and restart your shell after running 'conda init'.
```

conda认为我没有初始化环境，我脚本是在bash环境下运行的，我也用`conda init bash`初始化过。因此，问题肯定不是出在这里。

利用关键词"conda activate in bash script"检索，我找到了两种解决方法。

方法1: 在脚本中多加一句

```bash
source $HOME/miniconda/etc/profile.d/conda.sh
```

注意，我的conda是安装在家目录下的miniconda目录中，对于非家目录的安装方式，要修改 `$HOME/miniconda`。

方法2: 我们可以通过 `conda run` 来运行给定环境下的命令，假如，我们安装了一个环境rna-seq, 里面有一个程序叫做STAR, 我们可以随便写一个tmp.sh脚本，内容为

```bash
conda run -n rna-seq STAR --help
```

那么，此时运行 bash tmp.sh 就不会报错。也就是说，你并不是一定要用conda activate 启动环境，才能调用命令，你其实可以调用某个环境的给定指令。

方法2相对于方法1有个非常大的优势，那就是，如果你有多个不同python版本的环境，你不用担心写脚本的时候写了启动，但是忘了写退出。你只需要在原来的代码前加上一句， `conda run -n 环境名`。

- [calling conda source activate from bash script](https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script)
- [activating conda environment in within a shell script](https://askubuntu.com/questions/1218048/activating-conda-environment-in-within-a-shell-script)



