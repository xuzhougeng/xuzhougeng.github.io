---
title: C语言实战课-klib库学习
date: 2020-01-24 14:45:09.726
updated: 2020-01-24 16:18:00.166
url: /archives/c-practice-learn-klib
categories: 数据结构与算法
tags: C/C++
---

对于Python和R而言，使用一个已有的工具，通常都是先安装，然后看帮助文档学会怎么调用函数。那么对于C语言来说，我应该如何使用一个已有的轮子呢？

前些日子，我翻译了李恒大神一篇关于[seqtk](https://github.com/lh3/seqtk)代码介绍的博客，让我决定以seqtk为切入点，介绍下如何用别人造好的轮子。seqtk主要用到了两个头文件, khash.h和kseq.h, 同时我发现这两个文件在[klib](https://github.com/attractivechaos/klib)中，而klib提供了较为详细的介绍文档，因此我就直接去学klib了。

>关于klib名字中的k, 我想到了之前英语老师说，在英语中k和c的发音类似，例如kindle读起来就像candle, 因此这里的k以及函数名，结构体名里的k，都可以认为是c，表示这是一个c语言库。

## 第一步：阅读文档

一个优秀的项目必然会有一个优秀的文档对这个项目进行介绍，否则就只能自娱自乐，不被广泛使用。因此，学习klib的第一步就是阅读它的文档。我为了让自己更好的理解，因此通过翻译的方式让自己能够认真阅读。

为了实现范型的容器，klib广泛的使用了C宏。为了使用这些数据结构，我们需要先扩展出宏中的实例方法。于是，源代码就显得有点难以阅读以及不利于后续调试。由于C语言没有模版特性，实现高效范型编程也只能用宏。只有在宏的帮助下，我们才能写出范型容器，它在初始化之后才能在效率上和类型特异的容器竞争。其他一些库，例如[Glib](http://en.wikipedia.org/wiki/GLib), 使用`void *`类型实现容器。但是这种实现方法效率比较低，并且比klib使用更多内存。（参考这个[评测结果](http://attractivechaos.wordpress.com/2008/10/07/another-look-at-my-old-benchmark/)）

为了有效的使用klib，我们需要了解下它是如何实现范型编程。这里哈希表库（哈希表是一种键值对数据库，查询效率高）为例

```c
#include "khash.h"
KHASH_MAP_INIT_INT(32, char) //利用宏构建函数和数据结构
int main() {
	int ret, is_missing;
	khiter_t k;
	khash_t(32) *h = kh_init(32); //初始化
	k = kh_put(32, h, 5, &ret);   //插入key
	kh_value(h, k) = 10;   //根据key设置值
	k = kh_get(32, h, 10); //根据值找key
	is_missing = (k == kh_end(h));
	k = kh_get(32, h, 5);
	kh_del(32, h, k); //删除key
	for (k = kh_begin(h); k != kh_end(h); ++k) //遍历哈希表
		if (kh_exist(h, k)) kh_value(h, k) = 1;  //判断桶是否有数据
	kh_destroy(32, h); //删除哈希表
	return 0;
}
```

在这里例子中，第二行实例化了一个键类型为`unsigned`，值类型为`char`的哈希表，其中`m32`是哈希表的名字。所有和该名字相关的函数和数据类型都是宏，之后会进行解释。宏`kh_init()`初始化哈希表，而`kh_destroy()`则负责释放内存。`kh_put`插入键，并返回一个哈希表的迭代器(或者说位置)。`kh_get()`和`kh_del()`分别是获取值对应的键以及删除里面的元素。宏`kh_exist()`检测迭代器（或者说位置）是否在数据中。

看完代码之后，你会发现代码看起来不太像是有效的C程序，例如它缺少分号, 对一个函数赋值，`m32`没有预先定义。为了理解为什么代码是正确都，让我们更进一步都看下`khash.h`都源代码，代码骨架如下

```c
#define KHASH_INIT(name, SCOPE, key_t, val_t, is_map, _hashf, _hasheq) \
  typedef struct { \
    int n_buckets, size, n_occupied, upper_bound; \
    unsigned *flags; \
    key_t *keys; \
    val_t *vals; \
  } kh_##name##_t; \
  SCOPE inline kh_##name##_t *init_##name() { \
    return (kh_##name##_t*)calloc(1, sizeof(kh_##name##_t)); \
  } \
  SCOPE inline int get_##name(kh_##name##_t *h, key_t k) \
  ... \
  SCOPE inline void destroy_##name(kh_##name##_t *h) { \
    if (h) { \
      free(h->keys); free(h->flags); free(h->vals); free(h); \
    } \
  }

#define _int_hf(key) (unsigned)(key)
#define _int_heq(a, b) (a == b)
#define khash_t(name) kh_##name##_t
#define kh_value(h, k) ((h)->vals[k])
#define kh_begin(h, k) 0
#define kh_end(h) ((h)->n_buckets)
#define kh_init(name) init_##name()
#define kh_get(name, h, k) get_##name(h, k)
#define kh_destroy(name, h) destroy_##name(h)
...
#define KHASH_MAP_INIT_INT(name, val_t) \
	KHASH_INIT(name, static, unsigned, val_t, is_map, _int_hf, _int_heq)
```

`KHASH_INIT`是一个巨大都宏，他定义了所有都结构体方法。当这个宏被调用时，所有代码会通过C预处理器插入到程序中被调用都位置。如果一个宏被多次调用，那么代码会在程序中出现多个拷贝。为了避免哈希表中因为不同键-值类型导致都命名冲突，这个库使用[token concatenation](http://en.wikipedia.org/wiki/C_preprocessor#Token_concatenation), 也就是预处理器能够使用宏里的参数替换符号中部分内容。最后C预处理器会产生下面的代码，然后传递给编译器。（由于`kh_exist(h,k)`比较复杂，这里就不展开介绍）

```c
typedef struct {
  int n_buckets, size, n_occupied, upper_bound;
  unsigned *flags;
  unsigned *keys;
  char *vals;
} kh_m32_t;
static inline kh_m32_t *init_m32() {
  return (kh_m32_t*)calloc(1, sizeof(kh_m32_t));
}
static inline int get_m32(kh_m32_t *h, unsigned k)
...
static inline void destroy_m32(kh_m32_t *h) {
  if (h) {
    free(h->keys); free(h->flags); free(h->vals); free(h);
  }
}

int main() {
	int ret, is_missing;
	khint_t k;
	kh_m32_t *h = init_m32();
	k = put_m32(h, 5, &ret);
	if (!ret) del_m32(h, k);
	h->vals[k] = 10;
	k = get_m32(h, 10);
	is_missing = (k == h->n_buckets);
	k = get_m32(h, 5);
	del_m32(h, k);
	for (k = 0; k != h->n_buckets; ++k)
		if (kh_exist(h, k)) h->vals[k] = 1;
	destroy_m32(h);
	return 0;
}
```

从这个例子中，我们可以看到宏和C预处理器是klib中非常重要的部分。为什么klib那么快，一部分原因就是编译器在编译过程中已经知道键值对类型，因此他就能将代码优化到类特异代码对水平。一个用`void *`写的范型库无法达到这种性能。

在实例化过程中会有大量代码插入，这也提示我们为什么`C++`编译速度慢以及使用STL/boost编译的二进制文件很大。Klib由于代码量很小，并且各个组件相对独立，因此表现更加优异。插入上百行代码不会让编译变得很慢。

## 第二步：运行案例

介绍文档给了哈希表作为案例，我们就可以复制代码进行测试

```bash
git clone https://github.com/attractivechaos/klib
```

新建一个hash_test.c, 粘贴下列内容

```c
#include "klib/khash.h"
KHASH_MAP_INIT_INT(m32, char)        // instantiate structs and methods
int main() {
    int ret, is_missing;
    khint_t k;
    khash_t(m32) *h = kh_init(m32);  // allocate a hash table
    k = kh_put(m32, h, 5, &ret);     // insert a key to the hash table
    if (!ret) kh_del(m32, h, k);
    kh_value(h, k) = 10;             // set the value
    k = kh_get(m32, h, 10);          // query the hash table
    is_missing = (k == kh_end(h));   // test if the key is present
    k = kh_get(m32, h, 5);
    kh_del(m32, h, k);               // remove a key-value pair
    for (k = kh_begin(h); k != kh_end(h); ++k)  // traverse
        if (kh_exist(h, k))          // test if a bucket contains data
			kh_value(h, k) = 1;
    kh_destroy(m32, h);              // deallocate the hash table
    return 0;
}
```

编译和运行

```c
gcc -o test hash_test.c
./test
```

这个代码在运行时不会有任何输出，下一步我们可以尝试改造代码，让代码能够输出信息。

## 第三步: 改造代码

学会写代码的最快路径就是多写代码。以前我学习的时候，总觉得看一遍就够了，没有必要真的写一遍。真正开始写代码的时候，发现自己啥都不会。因此，为了能够学会使用khash.h, 我们按照自己的想法对代码做一些改造，看看代码的表现是否符合预期。下面代码，我写了一个循环，用于设置键值对。

```c
#include <stdio.h>
#include "klib/khash.h"
KHASH_MAP_INIT_INT(m32, char)        // instantiate structs and methods
int main() {
    int ret, is_missing;
    khint_t k;
    khash_t(m32) *h = kh_init(m32);  // allocate a hash table

	int i = 0;
	for ( i = 0 ; i < 10 ; i++ ){
        k = kh_put(m32, h, i, &ret);     // insert a key to the hash table
        if (!ret) kh_del(m32, h, k);
        kh_value(h, k) = i * i;             // set the value
	}

    for (k = kh_begin(h); k != kh_end(h); ++k)  // traverse
        if (kh_exist(h, k))          // test if a bucket contains data
			printf("%d\n", kh_value(h, k) );
    kh_destroy(m32, h);              // deallocate the hash table
    return 0;
}
```

假如你有个文件，也是键值对的关系，那么你可以读取这个文件，用文件里对内容对哈希表进行赋值。

## 第四步: 在项目中应用

学习klib的目的是为了解决我实际的问题，也就是写一个通用的fastq转换成fasta程序，这个程序能够自动解决行过长的问题。最终代码如下

```c
#include <stdio.h>
#include <zlib.h>
#include "klib/kseq.h"

KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{

        gzFile fp;
        kseq_t *seq;
        int l;
        if (argc == 1){
                fprintf(stderr, "Usage: %s <in.fasta|in.fasta.gz>\n", argv[0]);
                return 1;
        }

        fp = gzopen(argv[1], "r");
        seq = kseq_init(fp); // 分配内存给seq
        while( (l = kseq_read(seq)) >= 0){ //读取数据到seq中
                printf(">%s\n", seq->name.s); //
                printf("%s\n", seq->seq.s);
        }

        kseq_destroy(seq); //释放内存
        gzclose(fp);

        return 0;


}
```

这次的代码就特别的简洁，并且容易阅读。除了多了几行变量名定义和内存管理的代码外，其他时间就和写Python脚本一样，无非就是调用函数，`kseq_read`用于数据读取，`kseq_t`用于存放序列数据， 它的结构定义如下

```c
typedef struct __kstring_t {
        size_t l, m;
        char *s;
} kstring_t;
typedef struct {
        kstring_t name, comment, seq, qual;
        int last_char;
        kstream_t *f;
} kseq_t
```

这是我写的关于klib库的第一篇教程，往后还会不断用到这个库，用到的时候再写几篇进行介绍。

