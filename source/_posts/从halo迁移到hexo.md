---
title: 从halo迁移到hexo
date: 2023-08-18 15:16:49
url: /archives/from-halo-to-hexo
categories: 其他
tags: 个人博客
---

在距今很久的是2016年，我学会了购买域名，然后在一次课程中，我了解到了可以用github page搭建自己的网站，于是我倒腾了下Hexo，尝试这搭建了自己的博客。

但是由于水平不行，也不懂什么是github，最后也只是昙花一现，没能持久，所以，后来我用了简书、CSDN，将自己的技术博客发表在互联网上，让自己的内容更好被检索到。

然而，我内心还是有一个梦想的，我希望大家都可以通过访问包含我的名字拼音(xuzhougeng)的网站来看我的文章，后面机缘巧合下，我找到了一个博客搭建工具，叫做Halo。于是我就把自己的网站搬迁到Halo上，并坚持了多年。但是，由于我的服务器是在国外，就导致每次想把内容复制到微信公众号时，图片都会裂掉，我就需要花很多时间，做图片的重新复制工作。一种解决方案，就是我上个月重新买了一台国内服务器，然后重新部署了一个新版本的Halo。但是，让我难受的是，从1.x到2.x是一个不向下兼容的过程，就导致这个过程很痛苦。同时，每次登陆他的后台写东西也不是我中意的，我还是希望能够在本地写作，然后到处分发。

于是，我想了我的老朋友Hexo。我从之前的服务器上把数据导出成包含front-matter的markdown文件，然后将图片上传到腾讯的COS上，在把之前的图像URL做了替换，代码如下

```python
import re

cos_url = ""

os.makedirs('result', exist_ok=True)

for file in md_file:
    with open(file, 'r', encoding='utf-8') as f:
        content = f.read()
        # print(content)
        # print(re.findall(r'\(/upload/.*?\)', content))
        for img in re.findall(r'\(/upload/.*?\)', content):
            #print(img)
            # print(img[2:-1])
            content = content.replace(img, f'({cos_url}/' + img[2:-1] + ')')
            #print (content)
        # print(content)
        new_file = os.path.join('result', os.path.basename(file))
        with open(new_file, 'w', encoding='utf-8') as f:
            f.write(content)
```

参考，https://hexo.io/zh-cn/docs/github-pages.html 实现自动化部署，就得到了当前的xuzhougeng.top。

上面说到，为了保证图片的随处使用，我用了腾讯云的COS作为图床，用的是utools的图床插件，他的配置如下

![image-x0oiw2pb6n.png1692345113829.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/8/image-x0oiw2pb6n.png1692345113829.png)

但是目前的博客还是有点小问题，那就是这个tag不太好，同时内容也不是特别的全，缺少我简书里的内容，后续慢慢的优化吧。现在，让我们来发布这篇推文吧

```bash
hexo generate
git add *
git commit -m "从halo迁移到hexo"
# $Env:http_proxy="http://127.0.0.1:7890";$Env:https_proxy="http://127.0.0.1:7890"
# export http_proxy="http://127.0.0.1:7890";export https_proxy="http://127.0.0.1:7890"
git push -u origin main
```