---
title: MacOS原生Arm64架构的R语言环境配置-简明版
date: 2022-11-14 02:08:02.722
updated: 2022-11-14 02:08:02.722
url: /archives/macos原生arm64架构的r语言环境配置-简明版
categories: 
tags: 
---


在很早之前，我出过一期视频，介绍过如何在Arm64架构下（即M系列的芯片，包括M1，M2等）的Mac电脑上配置R语言分析环节。当时采用的是intel转译方案。只所以那样子做，是避免在环境配置上浪费不必要的时间，耽误了真正的学习过程。同时，我也写了[MacOS原生Arm64架构的R解决依赖编译问题](https://xuzhougeng.top/archives/how-to-compile-r-package-in-arm64-macos)给那些需要折腾的人。

不过，自从R进入了4.2, Bioconductor发布3.16之后，现阶段的M1/M2等用户就可以开始用原生的Arm版本的R了，因为无论是CRAN，还是Bioconductor都提供了预编译的arm版本的R包。
