---
title: 我的Rust第一课
date: 2021-08-27 14:01:50.324
updated: 2021-08-27 14:01:50.324
url: /archives/my-first-rust-class
categories: Rust
tags: 编程语言
---

在开始之前，先得做一些基础配置。

安装Rust工具链

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

我使用VScode进行开发，并安装了如下插件

- rust-analyzer
- rust syntax
- crates
- better toml
- rust test lens
- Tabnine

接着，尝试使用Rust构建第一个程序。我的目标是读取Fasta，统计其中的A, T, C, G各有多少个。

完成这个需求，有两个思路，第一个是调用成熟的库，第二个是自己造轮子写一个函数。

对于思路1，我使用关键词，Rust和Bio搜索到一个库叫做，[rust-bio](https://docs.rs/bio/0.37.1/bio/index.html)。 [https://rust-bio.github.io/](https://rust-bio.github.io/)

接下来，使用cargo new创建项目

```bash
cargo new fasta_count
```

第二步: 添加依赖。 编辑 Cargo.toml , 添加依赖信息

```Markup
[dependencies]
bio = "0.37.1"
```

第三步: 学习rust-bio文档，查找相关函数，然后在 `src/main.rs`中编写代码。

```rust
use bio::io::fasta;
use std::io;

fn main() {
    let  reader = fasta::Reader::new(io::stdin());
    let mut nb_a = 0;
    let mut nb_t = 0;
    let mut nb_c = 0;
    let mut nb_g = 0;
    for result in reader.records(){
        let record = result.expect("Error during fasta record parsing");
        for &base in record.seq() { //return sequence
            if base == b'a' || base == b'A' {
                nb_a += 1;
            } else if base == b'c' || base == b'C'{
                nb_c += 1;
            } else if base == b'g' || base == b'G'{ 
                nb_g += 1;
            } else if base == b't' || base == b'T'{ 
                nb_t += 1;
            }
        }
    }
    println!("A:{}, C:{}, G:{}, T:{}", nb_a, nb_c, nb_g, nb_t);
}


```

我使用标准库io::stdin从标准输入中读取数据，然后利用fasta::Reader的new出一个新的reader实例。对该实例的records对象进行遍历，得到fasta的中记录，record。最后对记录中序列进行遍历，得到不同碱基的数目。

第四步：编译代码。cargo会先将依赖环境下载到当前项目下，然后在进行编译。

```bash
cargo build 
```

编译结果存放在 target的debug目录下，名为fasta_count。我们运行程序，得到结果。

```Bash
cat TAIR10.fa | target/debug/fasta_count
# A:38223602, C:21551439, G:21528650, T:38177852
```

第一种思路的代码，是我看了极客时间「陈天·Rust编程第一课」中第一个示例"使用HTTP获取HTML然后保存成Markdown"后自己尝试的第一个和生信相关程序。对于第二种思路，我个人还是决定后续学习更多Rust基础再来尝试，因为我目前都不会用基础库来读取文件。

目前我的初体验是，Rust的语法更接近C/C++，它的编译器非常强大，能够很好的指出的我代码中的错误，便于检索debug。