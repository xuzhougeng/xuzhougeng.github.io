---
title: 在Windows中搭建一个基于WSL的Python数据科学环境
date: 2023-06-22 00:10:02.082
updated: 2023-06-22 00:18:29.794
url: /archives/building-a-wsl-based-python-data-science-environment-in-windows
categories: 数据科学
tags: 环境配置 | WSL
---


# 安装WSL，提供Linux环境

> 如果你发现后续的命令无法运行或者说软件商城中找不到，这可能意味着你的操作系统不符合要求。WSL安装要求 Windows 10 version 2004（Build 19041 ）及以上，或者是Windows11.

以管理员身份（也就是右击命令提示符）打开Windows下的CMD或PowerShell（后续，我们统一称之为**终端**）

![](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image_vySieCziap-9985fa1d46c14061828145444d4b3d33.png)

然后终端中，执行如下命令

```bash
wsl --install
```

中间**可能**会出现几次弹窗，需要用户进行确定，如下是命令运行时的界面。

![](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image_PQhwCz5ESA-7d76694de0dc4d52a5d0c42217ab783c.png)

或许？，你也可以通过**Microsoft Store**依序下载Windows Subsystem for Linux 和 Ubuntu（命令行的等价操作）。


![image_JZmsOMk0iJ](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image_JZmsOMk0iJ-03a42f94f0684bf8950bbaeff8ddf132.png)

最后安装完成后，需要你重启一下系统。

![](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image_odGFIUb3k--3c455a52782c41f88e01f729165a7d9e.png)

重启之后，会自动出现终端，让你输入用户名和密码，该用户还用sudo权限。

![image__lfZbSc_Fw](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image__lfZbSc_Fw-121185f65b994bc6a9498c9e9e9d7c0d.png)

如果你没有出现，或者自己手滑了，把它关闭了，这也问题不大，这意味着后续启动将以root用户启动。

PS：启动wsl的方式很简单，就是在终端里输入`wsl`并回车，或者在应用里搜索ubuntu

![image_vySieCziap](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image_vySieCziap-9985fa1d46c14061828145444d4b3d33.png)

如果你觉得权限太大，还是希望以普通用户的方式登录。我们可以在root下用`adduser` 增加一个用户【WSL环境内】

![image_HsYiTQ8y_i](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image_HsYiTQ8y_i-7888d4ea96cd4d1e9af5aefcd0bcfffa.png)


然后将刚才新建的用户添加到sudo组中，例如xzg

```bash
 usermod -a -G sudo xzg
```

最后在终端增设置我们刚才的用户为默认启动用户【直接打开的终端环境】

```bash
 ubuntu config --default-user xzg
```

后续启动wsl 的时候就不再是root，而是新建的用户，并且该用户也能使用sudo。

# 配置WSL环境

## 可选：系统代理

由于很多软件都需要从GitHub上下载，因此需要先为Linux系统配置一个系统代理。

先确定下自己的WSL版本

```bash
wsl -l -v
# 输出结果如下
  NAME      STATE           VERSION
* Ubuntu    Running         2

```

我的是WSL2，它是基于Hyper-V运行，因此与Windows系统在网络上是两台独立的机器。需要让WSL代理指向Windows的IP

```bash
# 这部分代码可以放在.bashrc中，启动wsl时自动配置环境
host_ip=$(cat /etc/resolv.conf |grep "nameserver" |cut -f 2 -d " ")
export ALL_PROXY="http://$host_ip:7890"
```

如果是WSL1，那么就只需要用下面这个代码

```bash
host_ip=127.0.0.1
export ALL_PROXY="http://$host_ip:7890"
```

之后，在Windows上打开Allow LAN，也就是允许本地区域网的访问。

![image_sSacJUuilI](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image_sSacJUuilI-718e4643ccda4335b8fea76c5490263a.png)

这样子后续数据下载走的都是代理。

## 可选：oh-my-zsh

之后，我们还需要简单的配置一个oh-my-zsh

```bash
# 安装zsh
sudo apt install zsh
# 安装oh-my-zsh
sh -c "$(curl -fsSL https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh)"

```

需要注意，我们后续安装conda之后，oh-my-zsh是没有conda配置的，因此我们还需要在bash环境中初始化下zsh相关的配置

```bash
conda init zsh
```

另外，配置在`.bashrc`里的一些内容，也应该需要自己手动加到`.zshrc`文件中

## 可选：NVIDIA驱动

因为我的电脑需要做深度学习相关的分析，所以需要配置深度学习的分析环境。

配置命令参考自：[https://developer.nvidia.com/cuda-downloads](https://developer.nvidia.com/cuda-downloads "https://developer.nvidia.com/cuda-downloads")

![image_yV0UklJJNt](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image_yV0UklJJNt-8056ceb7bd30461b8c542064d0e487ca.png)

只需要在WSL的Ubuntu中运行如下命令行

```bash
wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-wsl-ubuntu.pin
sudo mv cuda-wsl-ubuntu.pin /etc/apt/preferences.d/cuda-repository-pin-600
# 2.6G的数据下载
wget https://developer.download.nvidia.com/compute/cuda/12.1.1/local_installers/cuda-repo-wsl-ubuntu-12-1-local_12.1.1-1_amd64.deb
sudo dpkg -i cuda-repo-wsl-ubuntu-12-1-local_12.1.1-1_amd64.deb
sudo cp /var/cuda-repo-wsl-ubuntu-12-1-local/cuda-*-keyring.gpg /usr/share/keyrings/
sudo apt-get update
sudo apt-get -y install cuda

```

之后，只要在wsl的命令行中输入`nvidia-smi` 有输出就算成功。

```bash
nvidia-smi
# 我的输出如下
Wed Jun 21 23:26:37 2023
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 525.85.05    Driver Version: 528.24       CUDA Version: 12.0     |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|                               |                      |               MIG M. |
|===============================+======================+======================|
|   0  NVIDIA GeForce ...  On   | 00000000:01:00.0  On |                  Off |
|  0%   32C    P8    13W / 450W |   2884MiB / 24564MiB |      1%      Default |
|                               |                      |                  N/A |
+-------------------------------+----------------------+----------------------+

+-----------------------------------------------------------------------------+
| Processes:                                                                  |
|  GPU   GI   CI        PID   Type   Process name                  GPU Memory |
|        ID   ID                                                   Usage      |
|=============================================================================|
|    0   N/A  N/A        32      G   /Xwayland                       N/A      |
+-----------------------------------------------------------------------------+

```

如果系统上安装了pytorch，可以检查下，是否能够调用

```bash
import torch
torch.cuda.is_available()
# True 表示能够盗用
```

## 必选：安装conda

我们使用conda管理不同的环境，避免不同软件的依赖间冲突，导致环境崩溃。

打开你的wsl，我们下载mambaforge，使用mamba作为包管理器，并默认增加conda-forge源

```bash
# 下载mambaforge(需要不错的网络环境）
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
# 安装mamba
bash Mambaforge-$(uname)-$(uname -m).sh
# 安装时会有一些的选项，注意选择
```

之后，就可以用mamba隔离不同的环境。例如，我们需要安装一个pytorch环境（带GPU）

```bash
mamba create -n pytorch2.0.1 pytorch=2.0.1 torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia -y
# 大约要下载3.1GB
```

需要注意的是，mamba与镜像的兼容性不是特别好，因此不推荐使用mamba的时候修改镜像，而是使用官方源。

# 安装vscode

VScode从微软官方站点([https://code.visualstudio.com/](https://code.visualstudio.com/ "https://code.visualstudio.com/"))下载安装

![image_mH35pSo5-l](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image_mH35pSo5-l-ece8056d748e492ebe76cf2b8ed5e079.png)

## 配置插件

未安装插件的VSCode只是一款普通的编辑器，为了能让它满足我们后续分析数据的需求，我们需要为其配置一些插件。

我们需要先配置“Remote Development”，有了它，我们才能够连接到WSL系统。安装方式如下（其他插件也类似）

![image_VoS35suINB](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image_VoS35suINB-1725cb88fe0745f396cd68463cd86d81.png)

之后用快捷键shift + ctrl +p 调出，然后输入wsl，我们连接到WSL。

![image_VjdKtZP-0C](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image_VjdKtZP-0C-52f5fbda366b40648ec478a33afd8bda.png)

需要等待一会的后台初始化

![image_Iw04LJcB3A](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image_Iw04LJcB3A-8a50c36983e04b5480cb54d0639a627b.png)

## VSCode实现Jupyter

Jupyter 是一个开源的数据科学和机器学习工具，支持多种编程语言和平台。它允许用户在单个交互式笔记本中创建和共享文档，并且支持实时的代码编辑和运行。此外，Jupyter 还提供了一系列的库和扩展，可以用来完成各种数据分析、机器学习和科学计算任务。&#x20;

通过VSCode中上安装一些插件，也可以通过终端访问 Jupyter Notebook

-   Python: 使得VSCode能选择不同的Python内核
-   Pytlance：类型检查等功能
-   Jupyter: 在VSCode提供类Jupyter的界面

> 需要注意的是WSL系统下的插件和Windows下插件并不互通，要单独安装。

之后新建一个以.ipynb结尾的文件，就可以切换不同的Python内核

![image_t8wThT-sZh](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image_t8wThT-sZh-0f05902b41054071b963deddd2912e0d.png)

第一次运行代码时，可能会弹出如下提示，选择“安装”即可。

![image_RUBxu4j7TX](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image_RUBxu4j7TX-fa812a4b95b74def899479b821e4d770.png)

后续就跟网页版的Jupyter Notebook体验一致，你可以在其中画图，处理数据。

![image_ml67PloBt2](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/06/image_ml67PloBt2-3d430c2e05aa44759a22d219e100c99f.png)

那么以上，我们就算实现Windows系统下配置WSL，并使用VSCode作为我们的分析平台，用于处理数据

参考资料：

-   WSL的安装： [https://learn.microsoft.com/en-us/windows/wsl/install](https://learn.microsoft.com/en-us/windows/wsl/install "https://learn.microsoft.com/en-us/windows/wsl/install")
-   用户配置：[https://aka.ms/wslusers](https://aka.ms/wslusers "https://aka.ms/wslusers")&#x20;
-   WSL和VSCode的联动
-   WSL的系统配置：[https://learn.microsoft.com/zh-cn/windows/wsl/tutorials/wsl-vscode](https://learn.microsoft.com/zh-cn/windows/wsl/tutorials/wsl-vscode "https://learn.microsoft.com/zh-cn/windows/wsl/tutorials/wsl-vscode")
