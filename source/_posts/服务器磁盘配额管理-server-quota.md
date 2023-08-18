---
title: 服务器磁盘配额管理
date: 2019-12-27 10:06:48.616
updated: 2019-12-27 10:06:48.616
url: /archives/server-quota
categories: Linux
tags: 服务器
---


## 硬盘配额

Quota在是针对一个文件系统进行限制，而不是针对某个具体目录。也就是如果设置了`/dev/sda1`的quota，那么一个用户在这个`/dev/sda1`使用的数据量就会受限。如果你还有另外的`/dev/sda2`, 用户可以在超过`/dev/sda1`后继续使用`/dev/sda2`。

Quota针对整个文件系统的**限制项目**分为如下几个部分

- 容量限制或文件数目限制(block或inode)
- soft/hard: 超过soft时，用户会被警告，但是依然能操作，但是一旦超过hard，那么直接无法操作硬盘
- 宽限时间: 超过soft后，给用户一定时间进行数据转移

对于比较新的系统（例如CentOS7），建议在创建文件系统的时候就使用xfs，采用最新的quota策略。

### Ext4文件系统的Quota配置

第一步。增加文件系统支持

```bash
mount | grep /data2
# dev/sdb2 on /data2 type ext4 (rw,relatime,seclabel,stripe=64,data=ordered)
mount -o remount,usrquota,grpquota /data2
# /dev/sdb2 on /data2 type ext4 (rw,relatime,seclabel,quota,usrquota,grpquota,stripe=64,data=ordered)
```

第二步: 扫描文件系统并新建Quota的配置文件

```bash
quotacheck -avug
# -a: 扫描所有在`/etc/mtab`内含有quota支持的文件系统
# -v: 显示扫描过程
# -u: 新建aquota.user，按照用户进行扫描记录
# -g: 新建aquota.group, 按照用户组进行扫描记录
```

这一步可能出现的提示信息和解决方法

如果出现下面提示，说明内核支持xfs, 建议使用xfs的quota管理。

```bash
quotacheck: Your kernel probably supports journaled quota but you are not using it. Consider switching to journaled quota to avoid running quotacheck after an unclean shutdown.
```

如果出现下面提示，说明该磁盘正在独写，使用`lsof/fuser`查找可能的进程并关闭，之后重新运行。

```bash
quotacheck: Cannot remount filesystem mounted on /data2 read-only so counted values might not be right.
Please stop all programs writing to filesystem or use -m flag to force checking.
```

下面的提示会在首次执行时，可以忽略

```bash
quotacheck: Scanning /dev/sdb1 [/data3] done
quotacheck: Cannot stat old user quota file /data3/aquota.user: No such file or directory. Usage will not be subtracted.
quotacheck: Cannot stat old group quota file /data3/aquota.group: No such file or directory. Usage will not be subtracted.
quotacheck: Cannot stat old user quota file /data3/aquota.user: No such file or directory. Usage will not be subtracted.
quotacheck: Cannot stat old group quota file /data3/aquota.group: No such file or directory. Usage will not be subtracted.
quotacheck: Checked 44 directories and 23 files
quotacheck: Old file not found.
quotacheck: Old file not found.
```

第三步： 启动/关闭quota

```bash
quotaon -vug /data3
quotaoff -ug /data3
```

第四步: 编辑账户/用户组的限值和宽限时间

```bash
edquota -u 用户名
```

编辑表分为7个字段，修改的是第三个和第四个用于限制使用数据量(单位为KB)，修改第六个和第七个限制使用的文件数目(不常用)，0表示无限制。

和quota相关的Linux命令如下:

- quotactl
- quotaon/quotaoff
- edquota
- quotacheck
- repquota
- warnquota
- setquota

## XFS文件系统的Quota配置

如果使用xfs作为文件系统，那么就不能或者说没有必要采用上面的方法。xfs文件系统下，quota是文件系统的元信息，使用journaling提供更高层次的一致性保证，因此在管理上就会存在差别。

1. `quotacheck`对XFS系统没有效果。在mount时打开quota后，XFS会自动运行quotacheck做配额核算（quota accounting）
1. 没有必要在XFS文件系统下创建quota文件。aquota.user, aquota.group
1. 配额核算和限制增强之间存在区别。配额核算**必须**在mount的时候打开，但是限制增强可以在配额核算打开之后随意关闭和开启
1. 建议使用state监控XFS quota系统

挂载的时候打开quota

```bash
umount /mnt/disk1
mount -o uquota,gquota /dev/sbd1 /mnt/disk1
chmod 1777 /mnt/disk1
```

为了保证开机运行，还得修改`/etc/fstab`的第四列的`defaults`成`defaults,uquota,gquota`

> ext4可以通过`fstransform`转换成xfs

直接输入`xfs_quota`会进入交互模式

- help: 查看帮助信息, 例如help quota
- quota: 查看配额信息
- print: 列出设备
- df/free: 查看-bir(block,inode,realtime)信息，显示方式为 `-hN`,  `-f`指定文件
- quit/q： 退出交互模式

使用`-x`会进入专家模式， 之后通过`-c "命令"`执行命令, 举几个例子

```bash
# 显示使用情况
xfs_quota -x -c 'report -ugibh' /mnt/disk1
# 限制用户的文件数
xfs_quota -x -c 'limit -u isoft=5 ihard=10 user1' /mnt/disk1
# 限制用户的数据量
xfs_quota -x -c 'limit -u bsoft=1m bhard=2m user1' /mnt/disk1
```

其他可用命令

- path: 显示设备
- timer: 宽限时间，默认7天
- enable/disable: 开启或关闭配额增强
- off: 永久性关闭