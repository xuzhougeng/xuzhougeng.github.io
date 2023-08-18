---
title: 如何在Visual Studio上配置OpenCV 4.5.2用于C++项目开发
date: 2021-05-26 23:16:07.675
updated: 2021-05-26 23:19:51.423
url: /archives/visual-studio-configuration-of-opencv-for-cpp-development
categories: 数据科学
tags: C/C++
---

这是我第一次在号称宇宙最强的IDE Visual Studio上调用非标准库编译`C++`程序，中间遇到了很多的报错，让我怀疑人生。 

首先，我们需要从[Sourceforge](http://sourceforge.net/projects/opencvlibrary/files/opencv-win/)下载OpenCV, 目前最新版是4.5.2

![image.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-51a11bc996ed4923ba22cb81eb4891fa.png)

下载的exe文件双击之后会出现解压界面

![解压中](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-db97fe58403f431081103684dc9af996.png)

解压缩后得到如下文件

![解压内容](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-0ed0562fe5fa410c90109dfa60ffaad6.png)

我们需要依次系统属性->高级->环境变量，找到并设置设置环境变量Path

![环境变量](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-f9cec893255d4403a6cfe7295b6d1c9e.png)

如果不在环境变量中设置opencv中的bin路径，会出现如下报错 【由于找不到opencv_world452d.dll, 无法继续执行代码】

![error1](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-d79dd7544b154cc289b8c35bbfcd2387.png)

在Visual Studio 2019中使用OpenCV构建项目的流程如下

第一步：新建一个C++的控制台项目。

![新建项目](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-654548fdb7414748b468f464b1d3a27a.png)

第二步：配置项目的项目名称，例如opencv-test

![image.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-ad4394ad531d4675bd3d08c0fe732dff.png)

创建项目之后，先将Debug 改为x64，因为后续构建的是x64的项目。

![修改Debug](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-5234740fcd2d469caf3f511a4e14d6fc.png)

第三步：为了能让代码正常运行，我们需要配置这个OpenCV项目的所需的include和library路径，以及依赖的lib文件。

在菜单栏的项目(P)中选择属性(P)

![项目属性](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-151bb1e5a6d84d9b9a489adf37b5e54b.png)

选择VC++目录(VC++ Directories)的包含目录(Include Directories)，点击编辑

![image.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-135ed28733ac450093f5aa87bd5d22b0.png)

添加opencv的include的路径(头文件）

![include](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-0ac77d7dd74f43e4bcb8980d9555c3c6.png)

同样的操作，也用于库目录(Library Directories)的设置

![library](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-d83597426001403087700e7263bafe27.png)

添加library路径（路径比较深）

![library路径](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-1e3da7a7e9c243bca11ae07bc87d326f.png)

检查下，刚刚修改的include和library

![检查状态](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-d7c389f081cc4ab6b68ef0ff0fc9c24b.png)


此外还需要增加一个包括所有模块的lib文件。选择链接器(linker)中的输入(input), 接着编辑附加依赖项

![linker(https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-75a058ac43214ee9bfb9fb18f884df75.png)

添加**opencv_world452d.lib**，该文件包括opencv里的所有模块（lib文件和opencv版本有关，可从lib目录中确认）。

![符加依赖库](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-0ab452f73b074436bea8273df0612707.png)

如果不设置该变量，会出现【无法解析外部符号】的报错

![error2](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/05/image-2796e7747e1a4d09b3cb2d520973b6cb.png)

如果写错成 opencv_world452d.dll，就会出现【 LNK1104  无法打开文件"opencv_world452d.dll"】报错（之所以写成dll，是因为之前出现的dll找不到报错，导致我以为增加的是这个文件）。

如下是测试代码

```c++
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>

using namespace cv;
using namespace std;

int main()
{
  Mat image = Mat::zeros(600, 600, CV_8UC3);
  putText(image, "xuzhougeng", Point(100,300), FONT_HERSHEY_COMPLEX, 2.0, Scalar(100,200,200 ), 5 );
  imshow("Display Window", image);
  waitKey(0);
  return 0;
}
```

能正常出图表示能正确在项目中调用OpenCV。

参考资料：

- [https://docs.opencv.org/master/d3/d52/tutorial_windows_install.html](https://docs.opencv.org/master/d3/d52/tutorial_windows_install.html)
- [https://medium.com/@subwaymatch/opencv-410-with-vs-2019-3d0bc0c81d96](https://medium.com/@subwaymatch/opencv-410-with-vs-2019-3d0bc0c81d96)

