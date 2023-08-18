---
title: macOS的configd占了我好多内存
date: 2022-02-27 02:03:28.2
updated: 2022-02-27 02:03:28.2
url: /archives/macos-configd-use-too-much-memory
categories: Linux
tags: MacOS
---

在我没有启动多少应用的时候，macOS已经显示它使用了22.09GB内存。其中App内存是15.81GB, 我并没有打开那么多App.

![内存状态1](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2022/02/image-acb7830d509c43618a863043d689420e.png)

这估计跟configd有关，因为configd占用了20.55G内存。那么configd是真的占用了内存，还是就是声明自己会用到那么多内存呢？

我尝试着调用了比较多的内存，直接用了29.3Gb内存

```R
cols <- 8189
rows <- 320127
mat1 <- matrix(data = 0, nrow=320127, ncol = 8189)
print(object.size(mat1), unit="GB")
# 19.5 Gb
mat2 <- matrix(data = 0L, nrow=320127, ncol = 8189)
print(object.size(mat2), unit="GB")
# 9.8 Gb

```

此时看内存，发现rsession-arm64内存到了29.55GB, 而configd还是20.55GB，**已使用内存**仅仅比之前多了不到4GB。

![内存状态2](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2022/02/image-4c6c151de883487486c74b189992c246.png)

在我操作的前后，内存压力始终是绿色，说明系统没啥压力，那看来，有压力的显然是我。

我显然是之前用Windows系统，被某些软件的绿色气泡给吓到了，当气泡变红时，我就压力变大，就得点这个球释放内存。

而Unix系统的思路是调度所有内存，所有内存都要尽可能去干活，不要空着。

通过查阅文档，我们也可以知道configd是macOS特有的系统配置的后台驻留程序(System Configuration Daemon).configd管理本地系统的各类配置，维持数据映射到目标位置和当前系统状态，当数据发生改变时给应用程序发送信息，同时以  loadable bundles 形式托管一组配置代理。

最后，既然内存没压力，就不要自己给自己找压力，想着释放内存了。