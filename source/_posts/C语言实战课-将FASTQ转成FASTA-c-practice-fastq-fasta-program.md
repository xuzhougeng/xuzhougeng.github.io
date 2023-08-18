---
title: C语言实战课-将FASTQ转成FASTA
date: 2020-01-20 13:30:39.982
updated: 2020-01-20 13:33:57.626
url: /archives/c-practice-fastq-fasta-program
categories: 数据结构与算法
tags: C/C++
---

我这次学习C语言的目的非常单纯，就是尝试将C语言应用在日常的分析任务重，解决实际问题。既然如此，那么第一课就不是打印"hello world!" 了，毕竟说了那么多次世界你好，依旧写不好代码。

我们第一课直接就来处理一个实际的小需求，读取FASTQ，将其转成一个FASTA。处理这个问题和把大象放进冰箱里一样，都是分为三步，读取数据，处理数据，输出数据。其中第一步和第三步都是和文件打交道，而第二步考验的是对算法，数据结构和内存等有关知识对理解。

## 文件读写背景知识

我们应该如何读取数据呢？如果R语言，我会用`readLines("input.fasta")`直接读取所有的数据。如果是我以前用的python, 代码会是下面的样子

```python
fi = open("input.fasta", "r")
for line in fi.readline():
    print(line)
fi.close()
```

乍看起来，R好像比Python写的代码要少很多。但其实Python也可以一行完成，`[line for line in open("input.fasta", "r"]`。将数据输出到外部也差不多，基本上都能一行命令搞定。

总之，当你发现自己的代码可以少写了，其实是有人帮你简化了。

回归到C语言，所有的读写操作其实就是简单地将字节逐个移入程序中，将字节逐个从程序中移出而已，类似于水**流**动的感觉。

我们要让程序读取数据和输出数据，就相当于要搞一个水管，把水引过来。下面代码就是用`FILE *`声明了fi和fo这两根水管

```c
FILE *fi;
FILE *fo;
```

然后我们还需要把水管接到水流上。下面代码的`fopen`就是将我们之前对两根水管分别接到了输入文件和输出文件上，参数中的"r"表示读取(read), 而"w"表示写出(write)。

```c
fi = fopen("input.txt", "r");
fo = fopen("output.txt", "w");

//确保文件能够正常打开
//否则退出
if ( fi == NULL ){
    perror(argv[1]);
    exit(EXIT_FAILURE);
}
```

水管接上了，那么如何打开水阀，让水不断地流入流出呢？一种是逐个字符操作，一种是逐行操作。前面的方法更加基础，因为逐行读取就是不断的读取一个个字符，当碰到换行符时，就把前面的字符合在一起，以**字符串**的形式传递，输出是以行为单位，每一行也是一个个字符的写出。

- 逐字符函数: fgetc, fputc
- 逐行函数: fgets, fputs

通常读取和写出都伴随着数据处理，因此这部分放到实际数据处理进一步介绍。

最后关闭水龙头，我们需要用到`fclose()`

```c
fclose(fi);
fclose(fo);
```

## 动手搞事情

有了一些文件读写的基本概念之后，我们就可以开始动手写代码了。

根据fastq和fasta的格式定义，fasta是一行以`>`开头的序列表示，之后跟着N条序列。而fastq则是标准的4行，第一行是以`@`开头的序列标示符，第二行是序列，第三行以`+`开头，第四行是第二行序列对应的质量信息。处理过程描述如下

- 读取第一行
- 将第一个字符的`@`替换成`>`
- 输出第一行
- 读取第二行，输出第二行
- 读取第三行，不输出
- 读取第四行，不输出

下面实际的代码

```c
#include <stdio.h>
#include <stdlib.h>

#define MAX_LINE_LENGTH 1000

int main(int argc, char *argv[])
{
	FILE *fi;
	FILE *fo;

	if ( argc == 2) {
	    fi = fopen(argv[1], "r");
	} else if (argc == 3){
	    fo = fopen(argv[2], "w");
	} else {
		exit(EXIT_FAILURE);
	}

	char line[MAX_LINE_LENGTH];
	int count=0;

	while ( fgets(line, MAX_LINE_LENGTH, fi) != NULL ){
		if (count == 0){
		    line[0] = '>';
		}
		if (count < 2){
			if ( argc == 2){
			    fprintf(stdout, "%s", line);
			} else{
			    fprintf(fo, "%s", line);
			}
		}
		if (++count > 3){
			count = 0;
		}
	}

	fclose(fi);
	if (argc == 3) fclose(fo);

	return 0;

}
```

关于代码的讲解，我录制了专门的视频，<https://www.bilibili.com/video/av84245773/>

具体使用方法为

```bash
./fq2fa input.fq output.fa
```

目前代码代码比较简陋，至少存在下面这些问题

- 不支持gz压缩
- 不支持管道输入
- 没有文件类型检验，不能提示正确的错误

根据这些需求，我们就有了学习的目标，也就是后面课程的内容。