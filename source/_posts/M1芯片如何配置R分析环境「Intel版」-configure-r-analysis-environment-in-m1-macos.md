---
title: M1芯片如何配置R分析环境「Intel版」
date: 2022-02-01 12:17:34.54
updated: 2022-02-01 12:19:20.317
url: /archives/configure-r-analysis-environment-in-m1-macos
categories: R
tags: MacOS
---

> 视频地址: https://www.bilibili.com/video/BV1F3411E7td/

各位好，我是洲更，欢迎收看本期教程。

在这一期视频中，我将会根据我的个人经验介绍如何在M1芯片的Mac上配置R语言的分析环境。

自从苹果在2021年10月份发布搭载M1 Pro或者是M1 Max芯片的14寸和16寸的MacBook Pro之后，目前在苹果官方网站就只能够购买M1芯片的MacBook Pro了。也就意味着，往后我们购买到的Mac都是ARM加购，而不是原来intel的x86架构。

我在付完钱，等待2个星期之后，终于也在2021年的11月的最后一天，到我手上这台搭载M1 Max芯片的16寸的MacBook Pro 。

于是，我终于有设备可以让我介绍，如何在M1芯片的电脑上配置R语言的分析环境了。

> 插句题外话：不得不说，苹果的供应链系统真的强，只要你愿意等，你就能用它标的价格购买，而不像我之前想买的PS5和英伟达的显卡，要么就是几万人抢几百台，要么就得要多花个几百几千从其他店铺买。

在这次教程中，我会介绍如何在M1芯片的Mac电脑上搭建R的分析环境，主要包括如下几个部分

- R语言的安装
- RStudio的安装
- R包安装所需的环境，
- R包运行环境的配置
- R包安装

> 尽管，我主要是以M1芯片的Mac介绍如何配置R语言环境，但由于目前许多软件都没有跟上，所以下载软件大部分还是为intel芯片开发的，需要通过rosetta2转译后才能在M1芯片的mac使用，所以这篇教程也适用于搭载intel芯片的mac电脑。

## R语言的安装

首先，我们需要安装最核心的组件，也就是R语言。

我们先打开浏览器，搜索R，找到R语言的官网,[http://r-project.org](http://r-project.org), 选择Download R，国外用户选择主站，国内用户推荐选择清华镜像源

![R-Project](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2022/02/image-e7059bac0bf4463082a13f288ca6cbf4.png)

然后点击   Download R for macOS， 此处，我们就会有两个选择， 一个是intel版本的R，一个是arm64版本，也就是专门给M1芯片开发的R。

![R-下载](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2022/02/image-b035de5eb5164a88ba0766fcab4c2eda.png)

那么到底要选择哪个呢？

因为我这个视频是面向初学者的，比较基础，所以我无脑推荐使用Intel 64-bit的R语言，牺牲一点性能去避免不必要的折腾，是非常值得的。

后续，等到M1芯片相关的R语言环境变好了，我再出一期视频介绍。

下载完之后，我们只需要不断的下一步，就可以了。

## RStudio安装

说完了R语言的安装，RStudio的安装就非常的容易了，因为他没有选择，目前就只有一个intel的版本。

同样通过搜索找到Rstudio的官方网站，然后找到他的下载方式，然后点击下载Mac版本的Rstudio。下载地址为[https://www.rstudio.com/products/rstudio/download/#download](https://www.rstudio.com/products/rstudio/download/#download) 

![Rstudio](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2022/02/image-f0f9a435cdce4fd28bbdffdd444b71e9.png)

下载完成后，我们双击打开，将其拖动到Applications中。

然后打开启动台，等待Rstudio出现，然后点击它打开它。

这里会出现一个弹窗，问我们是不是需要安装一个命令行开发工具，我们选择取消。后续我们会专门去安装这个工具。

我们简单的打一个 print("hello world")测试即可。

在安装完R和Rstudio后，我还想给大家推荐一个非常好用的工具， [Rswitch](https://rud.is/rswitch/)

MacOS版本的RStudio是不支持切换不同的版本的R的，安装了这个工具后，我们就能够实现类似于Windows那种不同版本R语言切换的功能，

然而，这个软件比较小众，百度上我们暂时搜索不到（不过微软的必应可以），不过，这里我们选择直接输入网址， ([https://rud.is/rswitch/](https://rud.is/rswitch/))。

我们点击下载，等待下载完成，最后将其拖动到应用程序中就算安装完成了。

我们在启动台里将其打开，就可以在屏幕的右上角找到它。点击它，就能够用它来切换不同版本的R，以及用它来启动Rstudio。

![Rswitch说明](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2022/02/image-2e79b6309ba74b7ca55d77d51be1e1a0.png)

### 编译环境

在编译环境上，我们需要安装2个工具， xcode和gfortran。

xcode是苹果提供的开发工具，里面包括编译C和C++代码所需的工具。安装过程很简单，

- 首先，我们打开启动台，
- 然后，在其他里面，找到终端，并打开
- 使用快捷键cmd + + 就可以放大终端
- 在终端里输入 `sudo xcode-select --install`，敲回车键
- 此时提示我们需要输入密码，注意，为了保证安全，所以在输入密码的时候，是看不到我们输入内容，我们输入结束后，敲回车就可以确定了。
- 接着就有一个窗口弹出来，我们点击同意。

安装时间取决于网络环境，所以你的时间可能和我的不一样。我们等待它完成。

安装成功之后，我们在终端里面输入如下命令，如果没有提示找不到，就说明安装成功了。

```bash
which clang
which clang++
```



除了xcode外，R包编译还会用到gfortran。跟xcode不同，gfortran的安装需要分两种情况，一种是为intel版本的R准备的版本，一种是为M1芯片版本的R准备的版本。

那么gfortran从哪里下载呢？我们先在浏览器中输入这个地址，[https://mac.r-project.org/tools/index.html](https://mac.r-project.org/tools/index.html)，打开之后，就可以看到下载方式了。

因为我们安装的是intel 64位的R，所以我们需要选择intel的版本下载就行。

如果你发现自己无法打开，可以尝试[https://download.csdn.net/download/u012110870/78708981](https://download.csdn.net/download/u012110870/78708981)

intel的版本安装就是不断的下一步。


### 运行环境

接下来，我们来准备R语言的运行环境。

在运行环境上，我们需要安装两个工具，一个是Java,一个是，XQuartz, 

有些包在安装的时候看起来没啥问题，但是在运行的时候却会报错，比如说rJava，究其原因，就是因为我们没有安装java环境。

```R
library(rJava)

The operation couldn’t be completed. Unable to locate a Java Runtime.
Please visit http://www.java.com for information on installing Java.

错误: package or namespace load failed for ‘rJava’:
 loadNamespace()里算'rJava'时.onLoad失败了，详细内容：
  调用: fun(libname, pkgname)
  错误: JVM could not be found
此外: Warning messages:
1: In system("/usr/libexec/java_home", intern = TRUE) :
  运行命令'/usr/libexec/java_home'的状态是1
2: In fun(libname, pkgname) :
  Cannot find JVM library 'NA/lib/server/libjvm.dylib'
Install Java and/or check JAVA_HOME (if in doubt, do NOT set it, it will be detected)
```

那如何安装Java呢？这时候我们又需要借助于清华镜像了。

我们打开浏览器，搜索清华镜像，然后打开，选择 AdoptOpenJDK，我们选择最新的17，选择JDK，我们装的是intel版本的R，因此点击x64，点击mac，点击这个tar.gz结尾的文件进行下载。

下载结束之后，我们需要安装我们之前下载的安装包。这一步需要用到终端，请跟着视频输入命令。

这里，我选择将文件拖动到终端中，来获取文件路径

```Bash
sudo mkdir -p /opt/jdk_x64
sudo tar xf OpenJDK17U-jdk_x64_mac_hotspot_17.0.1_12.tar -C /opt/jdk_x64

```

按照上述步骤装完Java只完成一半，另外一半，我们需要在RStudio里操作。

我们打开Rstudio，输入如下命令 `file.edit("~/.Rprofile")` 打开R的配置文件，然后属于如下内容。由于你下载的版本可能和我的不一样，所以在输入过程中，要不断的使用tab 补全功能，确保路径正确。

```R
 Sys.setenv("JAVA_HOME"="/opt/jdk_x64/jdk-17.0.1+12/Contents/Home")
```

重启Rstudio后，R语言就会根据我们的配置来确定Java所在目录下。

除了java外，XQuartz，是后续安装的R包rgl不可或缺的一部分。如果没有安装XQuartz，即便安装了rgl，加载也会出现如下报错。

```R
> library(rgl)
Error in dyn.load(dynlib <- getDynlib(dir)) : 
  无法载入共享目标对象‘/Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/library/rgl/libs/rgl.so’：:
  dlopen(/Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/library/rgl/libs/rgl.so, 0x0006):
   Library not loaded: /opt/X11/lib/libGLU.1.dylib
  Referenced from: /Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/library/rgl/libs/rgl.so
  Reason: tried: '/opt/X11/lib/libGLU.1.dylib' (no such file), '/usr/lib/libGLU.1.dylib' (no such file)
错误: package or namespace load failed for ‘export’:
 loadNamespace()里算'rgl'时.onLoad失败了，详细内容：
  调用: rgl.init(initValue, onlyNULL)
  错误: OpenGL is not available in this build
此外: Warning messages:
1:   Loading rgl's DLL failed. 
  This build of rgl depends on XQuartz, which failed to load.
 See the discussion in https://stackoverflow.com/a/66127391/2554330 
2: Trying without OpenGL...
```

我们打开浏览器，输入它的官方站点: [https://www.xquartz.org/](https://www.xquartz.org/) 

XQuartz目前是2.8.1版本，我们选择下载，双击打开，开始安装。安装结束之后，mac会注销一下。之后，我们打开启动台，打开其他，就能看到它了。

## R包安装

在完成环境配置之后，我们终于可以安装R包了。

我们打开Rstudio，为了提供安装速度，我建议国内的用户修改一下镜像，提高安装速度，推荐国内的用户，按照如下方式修改修改。

首先，我们打开R的配置文件，然后我们打开清华镜像，点击CRAN旁边的问号，然后复制这一行到 R的配置文件里；接着我们找到Bioconductor，复制这一行到R的配置文件里。注意用cmd + s进行保存。

```r
file.edit("~/.Rprofile")

options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")


```

之后，我们关掉配置文件，并重启Rstudio。

我们新建一个R脚本。

我们首先安装一个tidyverse，来测试CRAN的包的安装。

输入install.packages("tidyverse"), 然后用快捷键 cmd + 回车运行。

首先蹦跶出一堆红字，这是它依赖的其他R包。

然后提示 有一个R包，它虽然有二进制版本，但源代码是最新的，问我们是不是要编译最新的。因为我们已经配置好了编译环境，所以这里输入yes。

后续，又出现了大量的红字，虽然很红，但是不是报错，都是正常信息。此时出现了黑字，这些事编译信息。

等他安装结束之后，我们输入 `library(tidyverse)` 来确认是否安装成功。 虽然有很多信息，但是没有看到error就说明成功了

接着，我们安装xlsx 来测试java环境。 没有报错说明，java环境正常。

然后是 rgl；我这里提示R出现fatal error，这并不是我们的环境配置出问题了，而是因为我的操作是在虚拟机里面执行的，是虚拟机性能不行，在加载时出错了。

接着，我们安装， BiocManager 用来安装 Bioconductor上的包。安装成功后，我们安装DESeq2，来测试Bioconductor，同样它也有很多依赖。

安装结束之后，提醒我们是否需要升级包，我这里输入n 表示不想升级，你们可以试试输入y。

最后我们试试安装github上的包，我们以ComplexHeat为例。我们输入如下命令， 运行的时候提示devtools没有找到，这说明我们没有安装这个包。

我们先去安装这个包； 这里有提示没有 这个函数， 原因是我输入出错了。改正就行了。

安装devtools后，我们再次运行之前的命令，这次就没有问题了。 我们加载也是成功的。

参考资料

- [https://mac.r-project.org/](https://mac.r-project.org/)
- [https://mac.r-project.org/tools/index.html](https://mac.r-project.org/tools/index.html)
- [https://kaopubear.top/blog/2021-07-06-macos-arm64-R-tips/](https://kaopubear.top/blog/2021-07-06-macos-arm64-R-tips/)
- [https://code2care.org/q/install-native-java-jdk-jre-on-apple-silicon-m1-mac](https://code2care.org/q/install-native-java-jdk-jre-on-apple-silicon-m1-mac) 