---
title: Windows中Visual studio使用Qt5的中文乱码问题
date: 2021-06-19 11:47:26.754
updated: 2021-06-19 11:47:59.148
url: /archives/windows-visual-studio-qt-jibberish
categories: C++
tags: QT | MSVC
---

Qt的中文输入至少包括如下两类：

1. 代码编辑器中的中文
2. 通过文件选择框打开的中文路径


对于第一种情况，代码编辑器中的中文跟**编辑器的编码**有关，因此可能会出现在编辑器中能够显示的中文，但是到了Qt的UI界面中就出现了乱码，例如

```c++
    QString status = "echo 测试 && pause";
    ui.textEdit->setText(status);    

```

![文本框中出现乱码](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/06/image-855e0b5aeba9497eb5b5cc56ddaa8ab2.png)

原因: Visual Studio中对中文的默认编码和QT对中文的默认编码不一样。也就是说，我们输入了X编码方案的中文，但是在编译阶段，Qt用了Y编码方案进行解码，从而导致源码阶段的乱码。

解决方案，就是设置好编码。

```C++
//声明来源是local 8bit
QString status = QString::fromLocal8Bit("echo 测试 && pause");
ui.textEdit->setText(status);
//用unicode字面量
QString status = fromUtf16(u"echo 测试 && pause");
ui.textEdit->setText(status);
```

对于第二种情况，QString获取到的中文就是它默认的编码方式，因此通过文件选择框获取的文件路径，如果包含中文，它也能够在Qt的UI界面正常展示。

但是如果直接转成string，然后调用system就会出现乱码，如

```c++
//乱码
QString status = QString::fromUtf16(u"echo 测试 && pause");
string test = status.toStdString();
system(test.c_str());
```

![cmd中的乱码](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/06/image-8e931980dca649a48e02ae1c05285a0b.png)

即，QString转string的时候存在编码上的问题。

解决方案, 将QString转成Loal8Bit才行，而不是直接使用toStdString，之后在构建成string，转成const char*才能避免乱码。

```c++
QString status = QString::fromUtf16(u"echo 测试 && pause");
string test(status.toLocal8Bit());
system(test.c_str());
```

参考资料

- [https://cloud.tencent.com/developer/article/1464366](https://cloud.tencent.com/developer/article/1464366)
- [彻底解决Qt中文乱码以及汉字编码的问题(UTF-8/GBK)](https://blog.csdn.net/libaineu2004/article/details/19245205)
- [借Qt中文乱码谈谈Coding中的编码问题](https://durant35.github.io/2016/02/02/programPearls_Qt_%E5%80%9FQt%E4%B8%AD%E6%96%87%E4%B9%B1%E7%A0%81%E8%B0%88%E8%B0%88Coding%E4%B8%AD%E7%9A%84%E7%BC%96%E7%A0%81%E9%97%AE%E9%A2%98/)