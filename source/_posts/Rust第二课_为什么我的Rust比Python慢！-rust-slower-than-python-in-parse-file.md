---
title: Rust第二课:为什么我的Rust比Python慢！
date: 2021-10-12 14:30:55.605
updated: 2021-10-12 14:30:55.605
url: /archives/rust-slower-than-python-in-parse-file
categories: Rust
tags: 算法
---

在[我的Rust第一课](https://xuzhougeng.top/archives/my-first-rust-class), 我写了一个程序对fasta中的ATCG进行计数。后面，我就想到一个非常常见的需求，对文件进行读取，统计行数，类似于 `wc -l`

下面是我写的第一个版本的代码, 我命名为myRead.rs

```rust
use std::io::BufReader;
use std::fs::File;
use std::env;
use std::io::BufRead;

fn main() -> std::io::Result<()> {

    let args: Vec<String> = env::args().collect();
    let filename = &args[1];

    let f = File::open(filename)?;
    let reader = BufReader::new(f);
    let mut line_num = 0;

    for _line in reader.lines() {
        line_num += 1

    }
    println!("{}", line_num);
    Ok(())
}

```

然后用rustc进行编译

```bash
rustc myRead.rs
```

接着我用一个记录cdna的fasta文件（约300MB)，进行测试，耗时约5.6秒

```bash
time ./myRead Homo_sapiens.GRCh38.cdna.all.fa
# 5.50s user 0.13s system 99% cpu 5.647 total
```


同样，我写了5行python脚本进行比较

```python
import sys
count = 0
for line in open(sys.argv[1]):
    count += 1
print(f"line number {count}")

```

python代码不到1秒就完成了任务

```bash
 time python ./read_file.py Homo_sapiens.GRCh38.cdna.all.fa
# 0.75s user 0.16s system 99% cpu 0.904 total
```

看到这个结果我直接震惊. Rust的运行速度居然比Python慢了6倍左右。经过高强度的检索，终于被我找到了靠谱的答案， [BufReader 100x slower than Python — am I doing something wrong?](https://users.rust-lang.org/t/bufreader-100x-slower-than-python-am-i-doing-something-wrong/1426)

总结下原因就是

1. 默认的优化不行，对于`rustc` 需要设置 `-C opt-level=2` 或者等价的 `-O`, 对于cargo则是设置 `--release`
2. `.lines()` 会为每一行都重新分配内存，因此不仅仅是处理UTF-8的问题。


根据第一个建议, 我重新用rustc编译了代码, 速度直接超过了Python

```bash
rustc -O ./myRead.rs
time ./myRead Homo_sapiens.GRCh38.cdna.all.fa
# 0.51s user 0.12s system 72% cpu 0.874 total
```

根据第二个建议，重新写了如下代码

```rust
use std::fs::File;
use std::env;
use std::io::Read;

fn main() {
    let args: Vec<String> = env::args().collect();
    let filename = &args[1];
    let mut file = File::open(filename).unwrap();
    let mut lines = 0;
    let mut buf = [0u8; 4096*32];
    while let Ok(num_bytes) = file.read(&mut buf) {
            if num_bytes == 0 { break; }
                lines += buf[..num_bytes].iter().filter(|&&byte| byte == b'\n').count();
    }
    println!("line number {}", lines);

}
```

使用 `rustc -O`编译后，运行速度又提升了2倍。

实际上，对于我这个Rust初学者，只需要记住Rust程序在编译的时候要**设置优化参数**, 进一步的优化代码，一时半会我还看不懂。






