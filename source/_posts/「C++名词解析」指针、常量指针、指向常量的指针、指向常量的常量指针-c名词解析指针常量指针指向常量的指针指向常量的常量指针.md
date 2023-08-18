---
title: 「C++名词解析」指针、常量指针、指向常量的指针、指向常量的常量指针
date: 2021-05-29 14:49:32.764
updated: 2021-05-29 23:28:04.768
url: /archives/c名词解析指针常量指针指向常量的指针指向常量的常量指针
categories: 
tags: 
---

指针(pointer), 一个有点特殊的变量，该变量记录着它所指向的对象的内存地址。

```c++
int a = 1;
int *ptr = &a;
```

常量指针(const pointer), 指针本质上是个变量，因此可以用 `const` 进行修饰，表示该指针记录的值是一个常量。常量指针在声明时需要初始化值。

```c++
int a=1;
int *const ptr = &a;
```

指向常量的指针(pointer refer to const)

```c++

int const a = 1;
int cont *ptr = &a;
```

指向常量的常量指针

```c++

int const *const ptr 

```