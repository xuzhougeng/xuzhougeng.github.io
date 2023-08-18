---
title: 「流畅的Python」关于函数参数，自省和闭包的一些记录
date: 2022-06-04 13:34:21.804
updated: 2022-06-04 13:34:21.804
url: /archives/流畅的python关于函数参数自省和闭包的一些记录
categories: 
tags: 
---

> 「流畅的Python」英文原版出版于2015年，2017年出版了中文版，我看到一些人推荐，于是就买了。但是遗憾的是那时候Python的道行太浅，发现每个字都认识，但是连起来就是读不懂。现在是2022年，由于疫情缘故，我终于有时间去深入学习一些Python知识。并且，经过几年的编程积累，我也终于有能力阅读这本书了。

“在Python中，函数是一等对象”。（最初看到这句话时，我脑海里想的是，编程语言还分个三六九等吗？）**一等对象**(first-class object, 又称第一类对象)是编程语言理论家定义的，符合下述要求的**程序实体**。

- 在运行时创建
- 能赋值给变量或数据结构的元素
- 能作为参数传给函数
- 能作为函数的返回结果

Python的函数之所以是一等对象，就是因为它满足上述条件，也就是能够实现下面这些操作

```python
def foo():
    print("foo")

#作为函数的参数
def bar(func):
    func()
    print("bar")

#作为函数的返回
def zoo(func):
    return func

f = foo #能赋值给变量
f()

ff = zoo(foo)
ff()
```

当一个函数能够接受函数作为参数，或者能够将函数作为结果进行返回，那么它就是一个**高阶函数**(higher-order function).因此上面代码的 `bar`, `zoo` 就是高阶函数。 

我们常用的sorted函数也是高阶函数, 我们可以通过传入自定义的函数来改变默认行为，下面的案例代码使得我们可以根据value对字典的key排序。

```python
d = {'a': 100, 'b': 1, 'c': 2}
# lambda是一个匿名函数, 这里是用来获取字典的值, 
sorted(d, key = lambda x : d[x] )
```

sorted中的key后面用到了lambda匿名函数，是一种非常简单的表达式，不支持try/catch, while,for等语句。由于lambda还会带来阅读上的困难，我们会更建议定义一个有名字的函数来替代。

Python函数支持非常灵活的参数处理机制，我们可以使用 `*`和`**`来展开可迭代对象，映射到单个参数中。例如

```python
#案例代码
def test_arguments(first, *middle, **last):
    print(f"first is: {first}") 
    
    for i,j in enumerate(middle):
        print(f"the {i} of middle is {j}")
        
    for k,v in last.items():
        print(f"key {k} is value {v}")

# middle捕获了b,c,d, 存成一个列表, last将k=1,k2=2存成一个字典, {"k":1,"k2":2}
test_arguments("a", "b", "c", "d", k=1, k2=2)

#first is: a
#the 0 of middle is b
#the 1 of middle is c
#the 2 of middle is d
#key k is value 1
#key k2 is value 2
```

因为Python一切皆对象，因此函数也是对象，既然是对象，我们就可以查看其属性（**函数自省**）。 Python使用 `dir`函数可以查看函数的所有属性, 如 `dir(sorted)`

> 这里的自省是一个计算机科学术语，是在程序运行时检查对象类型的一种能力，通常也称之为’运行时类型检查‘。我一开始还以为函数能自我反省，所以面对不理解的名词（包括后面出现的闭包），要果断去查，而不是由着自己乱想，。


