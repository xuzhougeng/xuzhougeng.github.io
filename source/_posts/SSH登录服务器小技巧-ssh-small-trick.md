---
title: SSH登录服务器小技巧
date: 2022-02-10 02:04:41.519
updated: 2022-02-10 02:05:51.707
url: /archives/ssh-small-trick
categories: Linux
tags: 小技巧 | SSH
---


最近切换到了MacOS平台进行办公，就不能用Windows下好用的XShell，用上了传统在命令行输入 `ssh -p port user@address`的方式进行登录了。

作为一个‘懒惰’的人，我肯定是要避免重复的运行登录命令了。回溯用过的命令进行复用是一种方式，但还是需要输入密码，所以我的操作方式如下

第一步: 通过编辑 `~/.ssh/config`文件， 为指定服务器增加别名

```conf
Host 别名
    HostName 服务器地址
    User 用户名
    Port 端口
```

这样子就能用 `ssh 目标服务器的别名`的方式登陆指定服务器，不必写后面的内容

第二步: 将本地的ssh公钥上传到目标服务器，实现免密登陆。

```bash
# 生成密钥
ssh-keygen
# 上传公钥
ssh-copy-id 目标服务器的别名
```

通过上述两步的配置，后续登陆服务器就就能节约不少时间，同时基于SSH的scp在服务器之间传送也不需要输入密码了。
