---
title: Rust第三课:模式匹配match语法
date: 2021-10-13 07:25:58.557
updated: 2021-10-15 07:20:39.708
url: /archives/rust第三课模式匹配match语法
categories: 
tags: 
---

> 我目前阅读Rust代码的一个不适应点就是Rust中的模式匹配match语法, 以及两个语法糖`if let`和`while let`，为了强化自己的记忆，于是便有了第三课。

所谓的模式匹配，就是对于不同情况采取不同的处理方法，

让我们想象一个应用场景，将一段DNA序列进行互补，即A->T, T->A, C->G, G->。 

对于Python而言，我们可以通过字典(dict)进行翻译

```python
dna = "ATCG"
trans_dict = {'A' : 'T', 'T': 'A', 'C':'G', 'G':'C'}
translated = ""
for base in dna:
    translated+=(trans_dict[base])

print(translated)
```

如果是C语言的话，我会用switch语法，如下

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(){
    char *dna  = "ATCG";
    int dna_len = strlen(dna);
    char *translated = (char*)malloc( sizeof(char) * dna_len);

    for (int i = 0; i < dna_len; i++){
        char base = dna[i];
	switch (base) {
            case 'A': 
                translated[i] = 'T';
                break;
            case 'T':
                translated[i] = 'A';
                break;
            case 'C':
                translated[i] = 'G';
                break;
            case 'G':
                translated[i] = 'C';
                break;
            default:
                translated[i] = 'N';
                break;
        }
    }
    printf("%s\n", translated);
    return 0;
}
```

Rust提供了match用于模式匹配。以一个碱基为例，我们输入'A', 希望翻译成'T', 下面是非常粗糙的代码

```rust
fn main(){
    let base = 'A';
    let trans = match base {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' =>  'C',
        _ => 'N',
    };
    println!("{}", trans);
```

其中match语法描述如下

```rust
match 变量名 {
    匹配类型1 => 代码1,
    匹配类型2 => 代码2,
    ....
    _         => 其他
}
```


match要求我们考虑所有的情况, 而一个utf-8字符肯定不只是ATCG这四个，因此其他情况需要用`_`指代，否则会报错。

对于一段序列, 我们可以写一个循环，直接在里面嵌套一个match语法。

```rust
fn main(){
    let dna = "ATCG".to_string();
    let mut translated = String::new();

    for i in dna.chars() {
        translated.push(
            match i {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                _   => 'N',
            }
            )
    }
    println!("{}", translated );

}

```

或者考虑将match匹配的部分封装到一个函数中

```rust
fn translate(base: char ) -> char {

    match base {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' =>  'C',
        _ => 'N',
    }
}

fn main(){

    let dna = "ATCG".to_string();
    let mut translated = String::new();

    for base in dna.chars() {
        translated.push(translate(base));

    }
    println!("{}", translated );
}
```

但是match不仅
