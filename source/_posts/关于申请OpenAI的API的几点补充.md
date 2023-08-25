---
title: 关于申请OpenAI的API的几点补充
date: 2023-08-25 08:41:57
update:
categories:
tags:
---

除了课程介绍如何使用第三方API外，可能大家也想申请一个OpenAI的官方账号。尽管在课程中，我们的演示过程非常顺利，中间没有出现意外，但考虑到不同学员的网络环境非常不同，所以这里做点补充。

## 网络环境

> 如果你使用[wildcard](https://bewildcard.com/i/ZHOUGENG)，可以使用他们提供的美国家庭的浏览器环境。

首先，一定要打开隐私模式或无痕模式或Edge浏览器的InPrivate模式。不同浏览器有不同的说法，以Edge为例

![InPrivate模式](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/8/image-dafflqhqfm.png1692924680509.png)

其次，打开全局代理。并且网络环境只能选择的非大陆地区，例如英美、新加坡等。部分线路由于使用人数过多，会被OpenAI或者谷歌邮箱拉黑，因此如果在操作过程中遇到这类问题，就需要切换线路。

注意，如果切换了线路，需要你重新打开一个隐私模式或无痕模式或Edge浏览器的InPrivate模式

关于网络环境的提供商，我个人使用的是PandaFan，它的价格比较贵，因此也就比较稳定。只不过，目前用的人好像多起来了，估计用于注册OpenAI和谷歌邮箱，可能不太好用，只能出现问题就换个地区了。

由于GFW的存在，因此PandaFan域名有很多个 https://pandafan.website/,https://pandafan.soccer/, https://pandafan.bid/ ,任选一个就行。（如果注册的时候想填写邀请码，可以考虑写185780）

## 谷歌邮箱申请

> 虽然OpenAI的账号不是所有人都需要，但是谷歌邮箱确实很实用。

如果自己有多年前申请的谷歌账号，那是最好的。如果没有，就需要单独申请一个。

部分同学可能在申请过程中，不需要借助短信验证。

如果需要短信验证，就需要使用到 https://sms-activate.org/cn。 

注册的谷歌账号，一定要注意，要先用一周（每天能够登录下），确保这个邮箱，不会被封号（俗称养号，不要被当做机器注册号），才能确保注册OpenAI后，能够一直使用谷歌邮箱登录。

## OpenAI注册

如果使用自己的环境在短信验证阶段、绑卡阶段出现被ban的情况，就不要在同一环境下继续尝试，建议换个环境，有群友换了，7、8次才成。
OpenAI注册时要用到短信验证，从下面中选择一个方法

方法1： https://sms-activate.org/cn 选择英、美地区，非虚拟的手机号>

![image-6k3guac8v9.png1692925772336.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/8/image-6k3guac8v9.png1692925772336.png)

方法2：使用bewildcard提供的5条免费短信（因为是免费，所以可能资源紧缺）。

## OpenAI绑卡付费

> 如果使用wildcard的环境也失败，联系下wildcard的客服进行解决。绑卡失败最大的问题就是网络环境，以及卡段问题。我当初就在这一步失败了很久。因此，课程中才介绍了使用第三方使用OpenAI的API方法。

在决定绑卡付费前，请先明白 https://chat.openai.com/ 的使用是不收钱的（当时要求网络环境），所以，如果不是刚需官方的api，可以先用一会 https://chat.openai.com/。

OpenAI的API的使用是收费的，并且由于常年被薅羊毛，新用户的免费额度一降再降，可能部分用户注册后连5美元都没有，甚至，你需要先存钱才能够调用，这都是正常的。因此，为了能够正常使用OpenAI的API付费，需要绑定信用卡。

早年间，我在折腾OpenAI付费的时候，经常遇到下面这个提示：

![你的卡被拒绝](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/8/image-n615s0f6hs.png1692925853007.png)

![your card is declined](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/8/image-zn3mr2gyjw.png1692925867483.png)

后来，我发现导致这个情况出现的的两个原因，就是虚拟信用卡的发卡商和网络环境。
经过不断的筛选，确定了课程中使用的[WildCard](https://bewildcard.com/i/ZHOUGENG)。这是因为我发现，只有wildcar都可以自动扣费，其他都会有问题。

另外请大家在绑卡的时候，最好使用他们提供的远程环境，他们的环境相比于我们自己的代理对OpenAI而言，绑卡更容易被通过。

这也是wildcard提供的解决方法 https://help.bewildcard.com/zh-CN/articles/8075545-%E7%BB%91%E5%AE%9A-openai-%E6%8F%90%E7%A4%BA-your-card-has-been-declined-%E6%80%8E%E4%B9%88%E5%8A%9E

![网络环境](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/8/image-v1gu0w8us2.png1692925985829.png)

OpenAI绑卡之后，会每个月自动根据你的这个月的使用量给你扣费（每月5号出账，确保卡上优余额），因此绑卡之后，后续原则上就用不上内部环境了。如果OpenAI扣款失败，请联系发卡商的客服进行解决。

最后，如果不打算继续使用OpenAI的API服务了，强烈建议取消付费计划，避免账号被盗导致扣费。

![取消](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/8/image-udqpizcgzv.png1692926013697.png)

## API使用

大部分基于OpenAI API开发的工具默认地址是OPENAI官方的API，即 https://api.openai.com/ 。但是，官方的API地址在使用时需要一定的网络条件。如果需要为亲朋好友、师长等提供帮助，建议用第三方，他们通常会为国内环境做优化。如果一定要使用官方的接口，可以参考 https://www.openai-proxy.com/ 自建转发服务器。

部分应用的API地址有特别的填写要求，例如utools的ChatGPT好友，需要加上 /v1/chat/completions。 请在使用的时候注意修改。
