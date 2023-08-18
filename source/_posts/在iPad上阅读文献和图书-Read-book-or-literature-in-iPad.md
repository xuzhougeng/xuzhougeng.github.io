---
title: 在iPad上阅读文献和图书
date: 2019-09-07 19:17:09.046
updated: 2019-12-17 18:11:13.117
url: /archives/Read-book-or-literature-in-iPad
categories: 其他
tags: Zotero | 文献 | 坚果云
---


长久以来，我用iPad读文献的操作一直是下面这几步

1. 用Zotero从网上抓取文献
1. 找出我的iPad连接线
1. 将iPad连接到电脑
1. 找到Zotero文献的所在目录
1. 打开iTunes的共享
1. 把PDF拖到iTnues的共享中
1. 打开iPad用GoodReader打开文献

每次插连接线都感觉自己太过原始。后来有了有了一部iPhoen手机，于是就变成了


1. 用Zotero从网上抓取文献
1. 打开电脑微信，把PDF发给自己
1. 打开iPhone的微信
1. 用iPhone的AirDrop功能将PDF共享给我的iPad
1. 打开iPad用GoodReader打开文献
1. 清理电脑微信和手机微信上的文献

手动清理是在是蠢。好在后来有了坚果云，就变成了

1. 用Zotero从网上抓取文献
1. 把文献保存到坚果云
1. 在iPad的GoodReader上同步文献

但是GoodReader这个软件只能作为一个阅读工具,它的PDF编辑功能太弱, 也不适合作为一个笔记工具

现在知道了PaperShip和MarginNote 3，于是就是

1. 用Zotero抓取文献
1. 用PaperShip同步文献，阅读文献
1. 假如文献信息量很大，用MarginNote 3做笔记
1. MarginNote 3做的笔记可以导出至印象笔记

~~不过PaperShip这个软件存在一些问题，会出现文件无法同步的问题，只能期待有更好的App了。~~感谢坚果云，目前PaperShip能够顺利进行同步了。

只要打开坚果云官网，登录坚果云账号，在zotero文件夹新建空白的lastsync.txt文件，**务必注意**的是，必须使用坚果云自带的新建文件工具来新建lastsync.txt文件，不能通过手动上传的方式。

![新建一个文件](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/12/image-396f9e8d19ab4096ab18603324e51a75.png)

接着，重新在Zotero客户端或者Papership验证服务器即可

![验证成功](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/12/image-cefd561f652241dfb287871b337e61a8.png)

> 如果我想看的是一本书，那么就是把书放到坚果云中，用MarginNote 3打开坚果云里的书做笔记。

## 利用ZotFile进行文献同步批注

还有一个解决方案，是ZotFile这款插件。我们只需要新建一个坚果云同步文件夹，Tablet，并在Zotfile的插件中加入该选项。

![同步文件夹](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/9/image-ca3d4cdc65654b88b891e9037b51f79f.png)

此时在Zotero处就会多出两个层次，存放送去批注的文件，以及批注完成的文件。

![Tablet Files](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/9/image-ac79c79640ac4af29c897195307e2a61.png)

我们将要阅读文献右击，在Manage Attachment 中选择 Send to Tablet

![送至平板](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/9/image-ac0822195ac34bdaabcfe836c7440d73.png)

因为我在GoodReader上设置坚果云中的Tablet作为同步服务器，因此就可以在GoodReader对文献进行批注。批注完成后，在Tablet Files(modified)中选择批注完成的文件，右击，在Manage Attachment 中选择 Get From Tablet

![取回修改后内容](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/9/image-2bfa2f008df64c7fab05e6965a83ccb1.png)

感觉这个方法也很不错，毕竟你回一次寝室看的文献也不会有那么多。

## 参考资料

- 思考问题的熊的博客: [Zotero入门学习路径](https://kaopubear.top/post/2019-09-11-howtolearnzotero/)
- 阳志平老师的博客: <https://www.yangzhiping.com/tech/zotero4.html>