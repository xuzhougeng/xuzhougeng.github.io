---
title: 如何让服务器不再龟速下载github的数据
date: 2023-09-02 18:57:44
update:
categories:
tags:
---

思路

1. 基于clash的代理，配置本地的代理转发
2. 配置环境变量，`http_proxy`和`https_proxy`


第一步：购买代理服务器，我推荐使用熊猫翻滚，即[PandaFan](https://PandaFan.soccer?r=185780/), 注册的时候可以考虑填写邀请ID, 185780。

第二步：安装clash，你可以考虑在服务器上直接下载，但是大概率会遇到龟速。百度网盘地址是，https://pan.baidu.com/s/1ERVR7_Y2hEilhqLKB4Gtrg?pwd=zgnh 提取码：zgnh .解压密码是 xuzhougeng.top。

```Bash
wget https://github.com/Dreamacro/clash/releases/download/v1.18.0/clash-linux-amd64-v1.18.0.gz
```

在服务器上的家目录中，解压缩clash，并赋予可执行权限。

```Bash
# 配置clash
mkdir clash
mv clash-linux-amd64-v1.18.0.gz clash
cd clash
gzip -d clash-linux-amd64-v1.18.0.gz
chmod +x clash-linux-amd64-v1.18.0
```

第三步，配置clash

clash的运行需要一个Country.mmdb文件，虽然会自动下载。但是速度比较慢，所以我也放在了百度盘中。按照如下方法进行手动配置

```bash
mkdir -p ~/.config/clash
mv ~/Country.mmdb ~/.config/clash
```

另外clash还需要一个yaml文件，配置代理规则，因此，你需要在clash的目录下新建一个config.yaml文件，内容如下

```yaml
port: 7890
socks-port: 7891
allow-lan: true
mode: Rule
log-level: info
external-controller: :9090
proxies:
- name: "trojan"
  type: trojan
  # interface-name: eth0
  # routing-mark: 1234
  server: 服务器地址
  port: 服务器端口
  password: 密码
proxy-groups:
  - name: science
    type: url-test
    proxies:
      - trojan
    url: http://www.gstatic.com/generate_204
    interval: 300
rules:
 - DOMAIN-SUFFIX,amazonaws.com,science
 - DOMAIN-KEYWORD,github,science
 - DOMAIN-SUFFIX,github.com,science
 - DOMAIN-SUFFIX,github.io,science
 - DOMAIN-SUFFIX,githubapp.com,science
 - DOMAIN-SUFFIX,githubassets.com,science
 - DOMAIN-SUFFIX,githubusercontent.com,science
 - DOMAIN-KEYWORD,github,science
 - DOMAIN-SUFFIX,github.com,science
 - DOMAIN-SUFFIX,github.io,science
 - DOMAIN-SUFFIX,githubusercontent.com,science
 - DOMAIN-SUFFIX,rawgithub.com,science
 - DOMAIN-SUFFIX,github.com,science
 - DOMAIN-SUFFIX,github.io,science
 - DOMAIN-SUFFIX,githubapp.com,science
 - DOMAIN-SUFFIX,githubassets.com,science
 - DOMAIN-SUFFIX,githubusercontent.com,science
 - GEOIP,CN,DIRECT
```

注意上面的server, port, password都需要自己配置。我们可以从pandafan的账号后台中获取

![配置方法](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/9/image-ay9p295qjj.png1693653031272.png)


第四步：在你需要从github上下载数据到时候， 我们先通过`./clash-linux-amd64-v1.18.0 -f config.yaml &`启动代理，然后配置环境变量`http_proxy`和`https_proxy`，如下

```Bash
export http_proxy=http://127.0.0.1:7890
export https_proxy=http://127.0.0.1:7890
```

那么后续从github上获取数据就非常方便了。例如从github上克隆singularity

```bash
git clone --recurse-submodules https://github.com/sylabs/singularity.git
```

速度如下

![git clone的速度](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/9/image-z4u68y52ab.png1693653372270.png)



参考资料

- [https://dreamacro.github.io/clash/](https://dreamacro.github.io/clash/)

