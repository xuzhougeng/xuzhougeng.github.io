---
title: 从typora转移到VNOTE
date: 2019-12-16 15:39:13.562
updated: 2019-12-16 15:39:13.562
url: /archives/typora-to-vnote
categories: 其他
tags: typora
---


因为之前一直用的是Typora进行写作，笔记管理的方式就是新建文件夹，也是下图这种情况。

![笔记本目录](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/12/image-254b1d8ed925480f9163a5c80f7cd0c1.png)

从轻度写作而言，Typora基本上够用了，但是它并不是一个笔记管理工具。熊的原话是，就本地文件的管理而言，Typora和VNote之间相差一到两个印象笔记。并且从LTF(List, Tag, Filter)角度来看，我只有list，而没有tag，更不要说filter了，因此查找笔记全靠记忆或者靠搜索引擎（我的学习笔记基本上都发布在互联网了）。

为了减少我转移笔记的手工操作，我分析了下VNOTE的笔记管理模式，发现核心文件就两个

- _v_imaes: 用于存放本地图片
- _vnote.json: 用于记录文件信息

由于并且我的图片存储习惯和VNOTE一致，也就是为不同笔记本建立单独的图片目录，我的是assets，而VNote是`_v_images`，因此我在添加笔记本的时候，只做了如下简单的更改。

![修改配置](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/12/image-6f2c53dfb35c4e31ab182b92590393de.png)

最终结果如下

![最终效果](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/12/image-d7b5eabbe57d4bc598fbc110d2b1df94.png)

因此，对于使用相对路径的我而言，我并不是从typora转移到VNote，而只是让我的之前的笔记本能够被VNote管理而已。
