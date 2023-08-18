---
title: 服务器硬盘管理笔记
date: 2019-09-05 22:50:20.564
updated: 2019-12-26 10:34:57.224
url: /archives/A-simple-disk-management-note
categories: Linux
tags: 服务器
---

# 服务器硬盘分区和挂载

## 文件系统简介

以ext4文件系统为例，设计的时候分为4个部分

- 超级块: 最开头部分，记录整个文件系统/分区中的文件数，`df`根据该信息进行统计磁盘占用
- 超级块副本: 对超级块进行备份，存在多个备份
- i节点(inode): 记录每个文件信息，名称/大小/编号/权限。 用`ls -i`获取i节点信息
- 数据块(datablock): 事实上的数据存储点，i节点以链接式结构记录数据块的信息，因此找到i节点，就能读取对应的数据。

由于`ls -l`获取的是i节点记录的数据使用的数据块个数，而`du`则是通过i节点获取实际大小, 所以`ls -l`和`du`显示的数据大小不同。

## RAID

RAID全称是Redundant Array of Independent Disks，也就是磁盘阵列，通过整合多块硬盘从而提升服务器数据的安全性，以及提高数据处理时的I/O性能。

RAID目前常用的是RAID5, 至少需要3块硬盘，其中一块硬盘用于奇偶校验，保证数据安全，其余硬盘同时读写，提高性能。此外，你还需要知道最原始的是RAID0，同时将数据读写到所有硬盘里，速度就变成了原来的N倍。RAID1至少需要两块盘，其中一块硬盘是另外硬盘的镜像。它不提高读写效率，只提高了数据安全性。RAID10是RAID0和RAID1的组合。

目前的服务器都配备了硬件RAID卡，因此在为服务器增加或更换硬盘时，需要**格外注意**，

- 插入的硬盘必须先做RAID才能被服务器识别
- RAID5中的所有硬盘应当被视为一个硬盘，不能只更换其中一部分

## 硬盘分区

fdisk只能对不多于2TB的硬盘进行分区

```bash
fdisk /dev/sdb
```

假如你的硬盘大于2TB，那么会输出如下信息

```bash
Welcome to fdisk (util-linux 2.23.2).

Changes will remain in memory only, until you decide to write them.
Be careful before using the write command.

Device does not contain a recognized partition table
Building a new DOS disklabel with disk identifier 0x405fb398.

WARNING: The size of this disk is 10.0 TB (10000294477824 bytes).
DOS partition table format can not be used on drives for volumes
larger than (2199023255040 bytes) for 512-byte sectors. Use parted(1) and GUID
partition table format (GPT).


The device presents a logical sector size that is smaller than
the physical sector size. Aligning to a physical sector (or optimal
I/O) size boundary is recommended, or performance may be impacted.

Command (m for help):
```

提示信息中的警告中，就建议"Use parted(1) and GUID  partition table format (GPT)."

因此，对于大于2TB的硬盘就需要用`parted`进行分区

```bash
parted /dev/sdb
```

输出信息如下

```bash
GNU Parted 3.1
Using /dev/sdb
Welcome to GNU Parted! Type 'help' to view a list of commands.
```

创建新的GPT标签，例如

```bash
mklabel gpt
```

设置单位

```bash
unit TB
```

创建分区, 比如我将原来的10T分成2TB和8TB

```bash
# mkpart PART-TYPE [FS-TYPE] START END
mkpart primary  ext4 0.00TB 2.00TB
mkpart primary  ext4 2.00TB 10.00TB
```

查看分区表

```bash
print
```

输出如下

```bash

Model: AVAGO MR9361-8i (scsi)
Disk /dev/sdb: 10.0TB
Sector size (logical/physical): 512B/4096B
Partition Table: gpt
Disk Flags: 

Number  Start   End     Size    File system  Name     Flags
 1      1049kB  2000GB  2000GB               primary
 2      2000GB  10.0TB  8000GB               primary

```

退出

```bash
quit
```

此时会提示"Information: You may need to update /etc/fstab." `/etc/fstab`用于设置开机硬盘自动挂载。如果硬盘被拔走了，而`/etc/fstab`没有修改，那么会就提示进行修复模式。

## 硬盘格式化

在挂载硬盘之前，需要先对磁盘进行格式化。使用的命令为`mkfs`, 使用`-t`指定文件系统，或者用`mkfs.xxx`，其中xxx就是对应的文件系统。文件系统有如下几类

- 传统文件系统: ext2, minix, msdos, fat, vfat
- 日志文件系统: ext3, xfs, ext4
- 网络文件系统: nfs

目前最流行的是ext4和xfs，足够稳定。其中xfs是CentOS7之后的默认文件系统。

```bash
mkfs.ext4 /dev/sdb1
mkfs.ext4 /dev/sdb2
```

## 硬盘挂载

之后用`mount`进行硬盘挂载，分别两种情况考虑

一种是新建一个文件路径，进行挂载。

```bash
mkdir data2
mount /dev/sdb2 /data2/
```

另一种是挂载一个已有目录，比如说临时文件目录`/tmp`挂载到新的设备中。

第一步: 新建一个挂载点，将原有数据移动到该目录下

```bash
mkdir /storage
mount /dev/sdb1  /storage
cp -pdr /tmp/* /storage
```

第二步: 删除原来的`/tmp`下内容

```bash
rm -rf /tmp/*
```

第三步: 重新挂载

```bash
umount /dev/sdb1
mount /dev/sdb1 /tmp
```

和mount相关的文件如下

- /etc/fstab: 文件系统表， 开机的时候会根据里面的记录进行硬盘挂载
- /etc/mtab: 记录着已经挂载的文件系统
- /etc/mtab~: 锁定文件
- /etc/mtab.tmp: 临时文件
- /etc/filesystems: 系统支持的文件系统

此外mount在挂载的时候还可以设置文件系统参数，例如是否支持磁盘配额，对应`-o`参数

## 总结

第零步: 检查服务器是否具备RAID阵列卡，如果有，则需要先为硬盘做RAID。

第一步: 使用`fdisk -l`检查硬盘是否能被系统检测到

第二步(可选):  假如需要**硬盘分区**，则用`fdisk/gdisk/parted`对硬盘划分磁盘

第三步: 使用mkfs进行磁盘**格式化**，有如下几种可选，

第四步: 用mkdir新建一个目录，然后用mount将格式化的硬盘挂载到指定目录下。卸载硬盘，则是`umout`

第五步: 修改`/etc/fstab `将硬盘在重启的时候自动挂载。**注意**: 如果硬盘不在了，则需要将对应行注释掉，否则会进入到emergency模式。
