---
title: R语言的字符串操作
date: 2019-07-28 10:25:03.264
updated: 2019-08-23 17:06:00.479
url: /archives/R语言的字符串操作
categories: R
tags: 正则表达式 | 字符串处理
---

R语言主要擅长于数值向量和矩阵操作，但是让他去做字符串操作也可以。

字符串的基本操作类型：

- 查找和替换
- 大小写转换
- 字符数统计
- 字符串连接和拆分

就我所知，有两套处理函数，一套是Hadley大神的stringr,一套是R自带的。

## stringr使用指南

stringr函数主要分为四类：

1. 字符操作：操作字符向量中的单个符 str_length, str_sub, str_dup

2. 添加，移除和操作空白符 str_pad, str_trim, str_wrap

3. 大小写转换处理 str_to_lower, str_to_upper, str_to_title

4. 模式匹配函数 str_detect, str_subset, str_count, str_locate, str_locate_all, str_match, str_match_all, str_replace, str_replace_all, str_split_fix, str_split, str_extract, str_extract_all 

### 单个字符的处理

字符长度`str_length`，等价于`nchar`

```
str_length("abc")
```

根据位置信息提取或替换字符, 类似于substr()

```
x <- c("abcdef","ghijk")
# 第三个
str_sub(x,3,3)
# 第二个到倒数第二个
str_sub(x,2, -2)
# 替换
str_sub(x,3,3) <- X
```

重复字符串，不同于rep

```
str_dup(x,c(2,3))
```

### 空白符

前后增加空白字符, str_pad()

```
x <- c("abc", "defghi")
str_pad(x, 10)
str_pad(x, 10, "both")
```

移除空白字符, str_trim()

```
x <- c("   b   ", "c   ", "   d")
str_trim(x)
str_trim(x,left)
```

更好的排版，让每一行的看起来一样长， str_wrap()

```
nature <- c("Nature Methods' Points of Significance column ",
"on statistics explains many key statistical and ","experimental design concepts. Other resources include an online",
" plotting tool and links to statistics guides from other publishers.")
cat(str_wrap(nature, width=40))
```

### 大小写转换

- 全部大写 str_to_upper, 类似于基础R的toupper()
- 全部小写 str_to_lower, 类似于基础R的tolower()
- title形式 str_to_title


### 模式匹配

功能： decect, locate, extract, match, replace, split

测试数据：

```
strings <- c(
    "apple",
    "219 733 8965",
    "329-293-8753",
    "Work: 579-499-7527; Home: 543.355.3679"
)

```

号码的正则形式:

```
phone <- "([2-9][0-9]{2})[- .]([0-9]{3})[- .]([0-9]{4})"
```

检测字符串是否符合特定模式， str_detect

```
str_detect(strings, phone)
```

提取字符串(全部内容）， str_subset

```
str_subset(strings, phone)
```

统计匹配次数，str_count

```
str_count(strings, phone)
```

定位首个匹配的位置,str_locate(返回matrix), 或所有匹配的位置, str_locate_all（返回list）

```
str_locate(strings, phone)
str_locate_all(strings, phone)
```

str_extract提取首个匹配，返回字符串向量, str_extract_all()提取所有匹配，返回list

```
str_extract(strings, phone)
str_extract_all(strings,phone)
str_extract_all(strings, phone, simplify=TRUE)
```

str_match提取首个匹配中()分组内的内容， str_match_all()则是全部，因此是list

```
str_match(strings, phone)
str_match_all(strings, phone)
```

str_replace()替代第一个匹配， str_replace_all替代所有匹配

```
str_replace(strings, phone, "XXX-XXX-XXXX")
str_replace_all(strings, phone, "XXX-XXX-XXXX")
```

str_split_fixed()返回固定数目，全部拆分，str_split（)可变拆分

```
str_split_fixed("a-b-c", "-")
str_spilt("a-b-c","-", n=2)
```

### 4类字符串描述引擎

1. 默认是正则表达式`vignette("regular-expression")

2. 逐byte固定匹配, fixed()

3. locale-sensive 字符匹配, coll()

4. 字符边界分析， boundary()

## 正则表达式练习

### 基本匹配

最简单的模式就是匹配某个完整的字符

```
x <- c("apple", "banana", "pear")
str_extract(x, "an")
```

如果需要忽略大小写, ignore_case =TRUE

```
bananas <- c("banana", "Banana", "BANANA")
str_dectect(bananas, regrex("banana", ignore_case=TRUE))
```

可以用点`.`匹配任意字符，但是默认不包括`\n`,需要用dotall=TURE开启

```
str_detect("\nX\n", ".X.")
str_detect("\nX\n", regex(".X.", dotall=TRUE))
```

### 转义（R中一坑）

正则表达式中有一些是特殊字符，比如说刚才的顿号，因此为了匹配这些特殊字符，我们需要对其转义。在Linux命令行里，转义用的是`\`, 所以可以直接用`\.`.

但是R里面的坑就出现了，我们用字符表示正则表达式，`\` 在字符里被用作转义符号。然后我们需要先把`\`转义成字符，然后才能进一步转义,`\\.`

如果要匹配`\` ,就需要`\\\\`，不可以思议，难以释怀，不知道被坑了多少次。

或者你用`\Q...\E` 类似于Python的 r'....'， 原意匹配

### 特殊字符

一些比较常用的字符匹配

- \d: 任意数字， \D：任意非数字.
- \s: 任意空白字符，\S：任意非空白字符
- \w: 匹配单词
- \b: 匹配字符边界， \B：非字符边界
- [abc], [a-z], [^abc], [\^\-]：匹配字母，和不匹配

在R里面，需要对"\"进行转义，所以上面的\在R里都要写成，\\
下面是一些预编译好的字符集，顾名思义

- [:punct:]
- [:alpha:]
- [:lower:]
- [:upper:]
- [:digit:]
- [:xdigit:]
- [:alnum:]
- [:cntrl:]
- [:graph:]
- [:print:]
- [:space:]
- [:blank:]

### 或

匹配abc或def

```
str_detect(c("abc","def","ghi"), "abc|def")
```

### 分组

匹配grey或gray

```
str_extract(c("grey","gray"), "gr(e|a)y")
```

分组可以用\1, \2进行提取，

### 定位

- `^`: 字符串开始， 如`^a`
- `$`: 字符串结束， 如`a$`

如果字符串有多行，那么就需要regex(multiline=TRUE)。此时，

- \A: 输入开头
- \z: 输入结尾
- \Z: 头尾

### 重复

- ?: 0或1
- +: 大于等于1
- \*： 大于等于0
- {n}: n次
- {n,m}： n到m次
- {n,}: 大于那次

默认是贪婪模式， 在上述字符后添加"?" 则为非贪婪模式。


PS: 下面是R语言自带的字符处理函数，我已经放弃他们了。

## 基础R包函数 

`nchar()`: 函数返回字符串长度
`paste()`, `paste0()`: 连接若干个字符串
`sprintf()`：格式化输出，下面举例

```
sprintf("%f", pi)
sprintf("%.3f", pi)
sprintf("%1.0f", pi)
sprintf("%5.1f", pi)
sprintf("%05.1f", pi)
sprintf("%+f", pi)
sprintf("% f", pi)
sprintf("%-10f", pi) # left justified
sprintf("%e", pi)
sprintf("%E", pi)
sprintf("%g", pi)
sprintf("%g",   1e6 * pi) # -> exponential
sprintf("%.9g", 1e6 * pi) # -> "fixed"
sprintf("%G", 1e-6 * pi)
```

`toupper()`: 大写转换
`tolower()`: 小写转换
`substr()`: 提取或替换一个字符串向量的子串

```
x <- "abcde"
substr(x,1,2)
# ab
substr(x,1,2) <- 2333
# 233cde
```

上面都是一些普通的函数，很好理解，下面都是一些和正则表达式相关的函数，如grep, grepl, regexpr, gregexpr, sub, gsub, strsplit
因此必须介绍一下R语言的正则表达式写法了。

- R语言是用的扩展正则表达式（Extended Regular Expressions)
- 元字符：\ | ( ) [ { ^ $ * + ?
- 非元字符转义后：\a as BEL, \e as ESC, \f as FF, \n as LF, \r as CR and \t as TAB
- 一些定义字符集合[:alnum:], [:alpha:], [:blank:], [:cntrl:], [:digit:], [:graph:], [:lower:], [:print:], [:punct:], [:space:], [:upper:],[:xdigit:]
- 找出“组”字符串 
- 默认是贪婪模式，可以通过用?改变为非贪婪模式

这些是基本知识，可以百度到每个字符的具体解释，或者看文档`?regexp`
不说基础知识了，看下应用吧。我常用的操作一般是找到某个字符串，或者对字符串进行替换

比如说，我想找到所有以P开头，且不是P结尾的字符，

```
test <- c("Python", "Perl", "PHP", "JAVA", "C", "C++")
grep("^P.*?[^P]$", test)
[1] 1 2
grep("^P.*?[^P]$", test,value=TRUE)
[1] "Python" "Perl"
grepl("^P.*?[^P]$", test)
[1]  TRUE  TRUE FALSE FALSE FALSE FALSE
regexpr("^P.*?[^P]$", test)
[1]  1  1 -1 -1 -1 -1
attr(,"match.length")
[1]  6  4 -1 -1 -1 -1
attr(,"useBytes")
> gregexpr("^P.*?[^P]$", test)
[[1]]
[1] 1
attr(,"match.length")
[1] 6
attr(,"useBytes")
[1] TRUE

[[2]]
[1] 1
attr(,"match.length")
[1] 4
attr(,"useBytes")
[1] TRUE
```

其中`grep()`默认是返回下标，如果设置value=TRUE，则返回字符串，`grepl()`返回是否配对的逻辑判断， `regexpr`则是返回匹配范围，如果不匹配结果是-1，`gregexpr`和前者功能一致，只不过返回的是列表形式。
**注**：忽略大小写ignore.case = TRUE

现在我想把`C++`替换成`C--`。我先试着找到`C++`

```
> grep("\+\+",test)
错误: 由""\+"开头的字符串中存在'\+'，但没有这种逸出号
```

什么情况，为什么`\+`不能把`+`这个元字符转义？难不成`+`在R里面不是元字符？我测试下

```
grep("++",test,value=TRUE)
Error in grep("++", test, value = TRUE) : 
  正规表现’++'不对，原因是'Invalid use of repetition operators'
```

啊！看来`+`还是元字符，难道是`\` 叛变革命了，我试试看。

```
> grep("\23","test\23",value=TRUE)
[1] "test\023"
> grep("\\23","test\23",value=TRUE)
Error in grep("\\23", "test\023", value = TRUE) : 
  正规表现’\23'不对，原因是'Invalid back reference'
```

看来`\` 是主要任务是把非元字符转义，如果想把元字符转义成普通字符，只能是`\\`元字符

```
grep("\\+\\+",test,value=TRUE)
[1] "C++"
```

回到我们之前的替换任务sub只对第一个匹配进行替换，gsub对所有匹配替换。

```
 sub("\\+\\+","--",test)
[1] "Python" "Perl"   "PHP"    "JAVA"   "C"      "C--" 
```

最后还可以用`strsplit`对字符串进行分割，返回的是一个列表

```
x <- c(as = "asfef", qu = "qwerty", "yuiop[", "b", "stuff.blah.yech")
strsplit(x, "e")
```


