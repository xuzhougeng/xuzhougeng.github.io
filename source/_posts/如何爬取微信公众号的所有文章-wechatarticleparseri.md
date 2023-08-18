---
title: 如何爬取微信公众号的所有文章
date: 2020-08-25 01:17:41.2
updated: 2020-08-25 01:18:37.893
url: /archives/wechatarticleparseri
categories: Python
tags: 爬虫
---


## 准备阶段

为了实现该爬虫我们需要用到如下工具

- Chrome浏览器
- Python 3 语法知识
- Python的Requests库

此外，这个爬取程序利用的是微信公众号后台编辑素材界面。原理是，当我们在插入超链接时，微信会调用专门的API（见下图），以获取指定公众号的文章列表。因此，我们还需要有一个公众号。

![image20200817090022286.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/08/image-20200817090022286-516f2abf9ee641638ed1315018eca56e.png)

## 正式开始

我们需要登录微信公众号，点击素材管理，点击新建图文消息，然后点击上方的超链接。

![image20200825110724858.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/08/image-20200825110724858-7684be606d6f4a629b037408487651b8.png)

接着，按F12，打开Chrome的开发者工具，选择Network


![image20200825110848472.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/08/image-20200825110848472-7946debc2b5a45bb8046e1c33d7d895f.png)

此时在之前的超链接界面中，点击「选择其他公众号」，输入你需要爬取的公众号（例如中国移动）

![image20200825111108639.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/08/image-20200825111108639-84afcf0c0ff74dc88aaa67f0c9b83dfa.png)

此时之前的Network就会刷新出一些链接，其中以"appmsg"开头的便是我们需要分析的内容

![image20200825111050522.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/08/image-20200825111050522-9152fdd915004685ae49ca784746e6a4.png)

我们解析请求的URL

```html
https://mp.weixin.qq.com/cgi-bin/appmsg?action=list_ex&begin=0&count=5&fakeid=MzI1MjU5MjMzNA==&type=9&query=&token=143406284&lang=zh_CN&f=json&ajax=1
```

它分为三个部分

- https://mp.weixin.qq.com/cgi-bin/appmsg: 请求的基础部分
- `?action=list_ex`: 常用于动态网站，实现不同的参数值而生成不同的页面或者返回不同的结果
- `&begin=0&count=5&fakeid`: 用于设置`?`里的参数，即begin=0, count=5

通过不断的浏览下一页，我们发现每次只有begin会发生变动，每次增加5，也就是count的值。

接着，我们通过Python来获取同样的资源，但直接运行如下代码是无法获取资源的

```python
import requests
url = "https://mp.weixin.qq.com/cgi-bin/appmsg?action=list_ex&begin=0&count=5&fakeid=MzI1MjU5MjMzNA==&type=9&query=&token=1957521839&lang=zh_CN&f=json&ajax=1"
requests.get(url).json() 
# {'base_resp': {'ret': 200003, 'err_msg': 'invalid session'}}
```

我们之所以能在浏览器上获取资源，是因为我们登录了微信公众号后端。而Python并没有我们的登录信息，所以请求是无效的。我们需要在requests中设置headers参数，在其中传入Cookie和User-Agent，来模拟登陆

由于每次头信息内容都会变动，因此我将这些内容放入在单独的文件中，即"wechat.yaml"，信息如下

```bash
cookie:  ua_id=wuzWM9FKE14...
user_agent: Mozilla/5.0...
```

之后只需要读取即可

```python
# 读取cookie和user_agent
import yaml
with open("wechat.yaml", "r") as file:
    file_data = file.read()
config = yaml.safe_load(file_data) 

headers = {
    "Cookie": config['cookie'],
    "User-Agent": config['user_agent'] 
}

requests.get(url, headers=headers, verify=False).json()
```

在返回的JSON中，我们就看到了每个文章的标题(title), 摘要(digest), 链接(link), 推送时间(update_time)和封面地址(cover)等信息。

> appmsgid是每一次推送的唯一标识符，aid则是每篇推文的唯一标识符。

![image20200817092802460.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/08/image-20200817092802460-0f93a4ed28ae4c8abbf72a56e360f395.png)

> 实际上，除了Cookie外，URL中的token参数也会用来限制爬虫，因此上述代码很有可能输出会是`{'base_resp': {'ret': 200040, 'err_msg': 'invalid csrf token'}}`

接着我们写一个循环，获取所有文章的JSON，并进行保存。

```python
import json
import requests
import time
import random

import yaml
with open("wechat.yaml", "r") as file:
    file_data = file.read()
config = yaml.safe_load(file_data) 

headers = {
    "Cookie": config['cookie'],
    "User-Agent": config['user_agent'] 
}

# 请求参数
url = "https://mp.weixin.qq.com/cgi-bin/appmsg"
begin = "0"
params = {
    "action": "list_ex",
    "begin": begin,
    "count": "5",
    "fakeid": config['fakeid'],
    "type": "9",
    "token": config['token'],
    "lang": "zh_CN",
    "f": "json",
    "ajax": "1"
}

# 存放结果
app_msg_list = []
# 在不知道公众号有多少文章的情况下，使用while语句
# 也方便重新运行时设置页数
i = 0
while True:
    begin = i * 5
    params["begin"] = str(begin)
    # 随机暂停几秒，避免过快的请求导致过快的被查到
    time.sleep(random.randint(1,10))
    resp = requests.get(url, headers=headers, params = params, verify=False)
    # 微信流量控制, 退出
    if resp.json()['base_resp']['ret'] == 200013:
        print("frequencey control, stop at {}".format(str(begin)))
        break
    
    # 如果返回的内容中为空则结束
    if len(resp.json()['app_msg_list']) == 0:
        print("all ariticle parsed")
        break
        
    app_msg_list.append(resp.json())
    # 翻页
    i += 1

```

在上面代码中，我将fakeid和token也存放在了"wechat.yaml"文件中，这是因为fakeid是每个公众号都特有的标识符，而token则会经常性变动，该信息既可以通过解析URL获取，也可以从开发者工具中查看

![image20200825123424376.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/08/image-20200825123424376-fb3e888613d14acf8f376e25a8fcf882.png)

在爬取一段时间后，就会遇到如下的问题

```python
{'base_resp': {'err_msg': 'freq control', 'ret': 200013}}
```

此时你在公众号后台尝试插入超链接时就能遇到如下这个提示

![image20200817102444104.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/08/image-20200817102444104-00823a33ae2e4acd989de9fc7d6e9d5d.png)

这是公众号的流量限制，通常需要等上30-60分钟才能继续。为了完美处理这个问题，你可能需要申请多个公众号，可能需要和微信公众号的登录系统斗智斗勇，或许还需要设置代理池。

但是我并不需要一个工业级别的爬虫，只想爬取自己公众号的信息，因此等个一小时，重新登录公众号，获取cookie和token，然后运行即可。我可不想用自己的兴趣挑战别人的饭碗。

最后将结果以JSON格式保存。

```python
# 保存结果为JSON
json_name = "mp_data_{}.json".format(str(begin))
with open(json_name, "w") as file:
    file.write(json.dumps(app_msg_list, indent=2, ensure_ascii=False))
```

或者提取文章标识符，标题，URL，发布时间这四列信息，保存成CSV。

```python
info_list = []
for msg in app_msg_list:
    if "app_msg_list" in msg:
        for item in msg["app_msg_list"]:
            info = '"{}","{}","{}","{}"'.format(str(item["aid"]), item['title'], item['link'], str(item['create_time']))
            info_list.append(info)
# save as csv
with open("app_msg_list.csv", "w") as file:
    file.writelines("\n".join(info_list))         
```

下一篇，将介绍如何根据每个文章的连接地址，来获取每篇文章的阅读量信息。

## 参考资料

- https://blog.csdn.net/kindred_joe/article/details/99289890
- https://blog.csdn.net/qq_28804275/article/details/82150874



最终代码如下，使用方法为`python wechat_parser.py wechat.yaml`

```python
import json
import requests
import time
import random
import os
import yaml
import sys

if len(sys.argv) < 2:
    print("too few arguments")
    sys.exit(1)

yaml_file = sys.argv[1]
if not os.path.exists(yaml_file):
    print("yaml_file is not exists")
    sys.exit(1)
    

with open(yaml_file, "r") as file:
    file_data = file.read()
config = yaml.safe_load(file_data)

headers = {
    "Cookie": config['cookie'],
    "User-Agent": config['user_agent'] 
}

# 请求参数
url = "https://mp.weixin.qq.com/cgi-bin/appmsg"
begin = "0"
params = {
    "action": "list_ex",
    "begin": begin,
    "count": "5",
    "fakeid": config['fakeid'],
    "type": "9",
    "token": config['token'],
    "lang": "zh_CN",
    "f": "json",
    "ajax": "1"
}

# 存放结果
if os.path.exists("mp_data.json"):
    with open("mp_data.json", "r") as file:
        app_msg_list = json.load(file)
else:
    app_msg_list = []
# 在不知道公众号有多少文章的情况下，使用while语句
# 也方便重新运行时设置页数
i = len(app_msg_list) // 5
while True:
    begin = i * 5
    params["begin"] = str(begin)
    # 随机暂停几秒，避免过快的请求导致过快的被查到
    time.sleep(random.randint(1,10))
    resp = requests.get(url, headers=headers, params = params, verify=False)
    # 微信流量控制, 退出
    if resp.json()['base_resp']['ret'] == 200013:
        print("frequencey control, stop at {}".format(str(begin)))
        break
    
    # 如果返回的内容中为空则结束
    if len(resp.json()['app_msg_list']) == 0:
        print("all ariticle parsed")
        break
        
    app_msg_list.append(resp.json())
    # 翻页
    i += 1

# 保存结果为JSON
json_name = "mp_data.json"
with open(json_name, "w") as file:
    file.write(json.dumps(app_msg_list, indent=2, ensure_ascii=False))
```