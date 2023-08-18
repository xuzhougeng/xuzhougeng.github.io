---
title: 一次路由器刷固件安插件的折腾之旅
date: 2021-08-08 23:04:16.12
updated: 2021-08-08 23:04:16.12
url: /archives/update-router-firmware-for-scientific-internet-access
categories: 其他
tags: 小技巧
---

故事的起因可以追溯到好几年之前，当时机缘巧合之下和[Leo@China](https://github.com/leoatchina)一起拼团买了一个二手斐讯路由器（关于这个路由器，还有一个金融理财事件在里面），第一次听说路由器刷固件的事情。

最近为了改善我们课题组的上网环境（人多设备多，旧的路由器容易断连），斥巨资购买了华硕路由器，型号为RT-AX86U，在给路由器设置时，发现界面有些眼熟，不由得让我想起来之前那台斐讯路由器的后台，于是我就想着，这台路由器是不是也能刷固件呢？

经过一波上网检索，我最终折腾成功，给这台路由器刷上了固件，离线安装了个别插件，如下是具体的教程。

> 如果没有特别的上网需求，不建议刷固件。

第一步，到KoolShare论坛的固件下载服务器， https://firmware.koolshare.cn/， 选择游客所能下载的最新固件。我购买的是2021年的RT-AX86U, 出场的固件版本就是386，因此我在 https://firmware.koolshare.cn/Koolshare_RMerl_New_Gen_386/RT-AX86U/ 里选择7月30日发布的386固件进行下载

![选择固件](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/08/image-8185fa73fbe34f7fa3145925d0907f23.png)

第二步：通过 router.asus.com 登录到路由器的后台，选择【系统管理】-> 【固件升级】， 上传下载的固件

![上传固件](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/08/image-53f126be31314376ba12e05650f3681a.png)

静静等待固件升级完成。

第三步：登录路由器，选择最后的 【软件中心】，到【未安装】中安装 shellinabox，启动shell命令行，运行如下命令

```bash
sed -i 's/\tdetect_package/\t# detect_package/g' /koolshare/scripts/ks_tar_install.sh
```

第四步：离线安装插件。访问 https://github.com/hq450/fancyss ，点击 fancyss_hnd ，下载其中的 .tar.gz文件。接着，打开路由器后台，选择最后的软件中心的离线安装，上传该 tar.gz 文件，等待安装完成。

第五步: 配置第四步里安装的插件，即可实现更加自由的网页访问。

最后，感谢KoolShare论坛和GitHub的hq450。 