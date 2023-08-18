---
title: C语言实战课-让程序能够处理gz文件
date: 2020-01-23 23:34:01.877
updated: 2020-04-27 22:19:27.688
url: /archives/c-practice-fastq-add-gz-support
categories: 数据结构与算法
tags: C/C++
---

之前的程序不能够处理压缩文件，而事实上为了节约空间，基本上fastq都会压缩成gz格式，因此这一课就是让程序能够支持压缩文件。

这里有两种思路，一种是利用管道，将之前的压缩文件通过`zcat`程序读取然后传递给我们的程序，另一种则是在程序中调用zlib库，让程序能够直接处理gz文件.

## 背景知识

**针对思路1**: 系统在每个C语言运行的时候都至少会提供三个流，标准输入(stdin)，标准输出(stdout)和标准错误(stderr).  在之前的程序我们就用到了stdout，用于将结果输出到屏幕上, 即`fprintf(stdout, "%s", line);`。 同样，我们可以修改`fgets(line, MAX_LINE_LENGTH, fi) `中的`fi`为`stdin`使得程序能够接受管道传递的数据。

**针对思路2**: 利用zlib读取gz文件并不复杂，只需要将原来的函数名前面加上或者改成gz。毕竟优秀的代码应该符合人的直觉，下面就是我们将要用到的几个函数

- 声明水管: gzFile
- 连接水管: `gzopen(const char *path, const char *mode)`
- 读取一行: `gzgets(gzFile file, char *buf, int len)`
- 输出到gz文件中: `gzprintf(gzFile file, const char *format, ...)`
- 关闭水管: `gzclose(gzFile file)`

我们以一个非常简单的代码进行展示，它会从gz文件中按行读取然后输出

```c
#include <stdio.h>
#include <zlib.h>

int main(int argc, char *argv[])
{

	gzFile fi;
	fi = gzopen(argv[1], "r");
	char string[1000];

	while ( gzgets(fi, string, 1000) != NULL) {
		printf("%s", string);
	}
	return 0;
}
```

它不是只能处理压缩文件，对于没有压缩的文件也是可以识别的，因此不需要写专门的逻辑判断语句，根据文件名来决定是否使用zlib相关的函数。

>  gzopen can be used to read a file which is not in gzip format; in this case gzread will directly read from the file without decompression.  When reading, this will be detected automatically by looking for the magic two-byte gzip header

由于zlib不是标准库，我们在编译程序的时候需要额外加上`-lz`参数，用于链接动态库。

## 动手搞事情

有了以上的背景知识后，我们就可以来修改我们之前的代码了.

zlib在大部分服务器上都属于默认安装，因此只需要多加一行`#include <zlib.h>`. 但是还是有一小部分服务器居然没有安装zlib，因此我们需要先下载编译（root权限可以直接用apt/yum进行安装）。

```bash
wget https://www.zlib.net/zlib-1.2.11.tar.gz
tar xf zlib-1.2.11.tar.gz
cd zlib-1.2.11 && ./configure && make -j 8
```

修改后的fq2fa.c的代码如下

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zlib.h"

#define MAX_LINE_LENGTH 1000

int main(int argc, char *argv[]){

	gzFile fi;
	gzFile fo;
	if (argc < 2 ){
		return -1;
	}
	fi = gzopen( argv[1], "rb");
	if ( argc == 3 ){
	    fo = gzopen (argv[2], "wb");
	}

	char line[MAX_LINE_LENGTH];
	int count = 0;
	while ( gzgets(fi,line, MAX_LINE_LENGTH) != NULL) {
		if ( count == 0) {
			line[0] = '>';
		}
		if ( count < 2){
			if ( argc == 2) {
				fprintf(stdout, "%s", line);
			}
			if ( argc ==3 ){
			    gzprintf(fo, "%s", line);
			}
		}
		if ( ++count > 3 ){
			count = 0;
		}

	}
	gzclose(fi);
	if (argc == 3) gzclose(fo);
	return 0;
}
```

编译方法为`gcc -lz -Lzlib -o fq2fa fq2fa.c`。关于代码的讲解和zlib的知识点可以看视频介绍, <https://www.bilibili.com/video/av84768292/>。

代码目前存在的问题是行缓冲是1000字符，对于三代测序数据，按照宣传可以有几M，一般也是几百k，因此当前程序处理三代程序肯定不行。一个方法是将`MAX_LINE_LENGTH`设置为更大的值，比如说设置为10M，或者是将其设置为一个参数，可以手动指定。

不过下一次，将会学习一个新的轮子[klib](https://github.com/attractivechaos/klib/)，用里面的函数来解决这个问题。