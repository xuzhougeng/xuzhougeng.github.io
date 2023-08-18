---
title: Centos的Slurm安装笔记
date: 2021-04-21 11:49:05.177
updated: 2021-04-21 11:49:05.177
url: /archives/install-slurm-in-centos
categories: Linux
tags: 小技巧
---

> 因为有一些软件必须要用Slurm，所以不得不在我的主机上配置slurm。

Slurm的安装依赖于root权限

### munge配置

```bash
wget https://github.com/dun/munge/releases/download/munge-0.5.14/munge-0.5.14.tar.xz
rpmbuild -tb --without verify munge-0.5.14.tar.xz
cd /root/rpmbuild/RPMS/x86_64
rpm -ivh munge-0.5.14-1.el7.x86_64.rpm \
    munge-devel-0.5.14-1.el7.x86_64.rpm munge-libs-0.5.14-1.el7.x86_64.rpm
```

创建密钥

```bash
sudo -u munge /usr/sbin/mungekey -v
# mungekey: Info: Created "/etc/munge/munge.key" with 1024-bit key 
```

生成的 munge.key 文件需要分发到所有的计算节点。

启动守护进程(daemon)

```bash
systemctl enable munge
systemctl start  munge
# 检查状态
systemctl status munge 
```

### 方法1: RPM安装

下载页面， [https://www.schedmd.com/downloads.php](https://www.schedmd.com/downloads.php)

因为是CentOS7, 因此我下载的是19.05版本。 而20.11可能不再支持Python2。

```bash
wget https://download.schedmd.com/slurm/slurm-19.05.8.tar.bz2
yum install pam-devel perl-Switch -y
rpmbuild -ta slurm-19.05.8.tar.bz2
cd /root/rpmbuild/RPMS/x86_64
rpm --install slurm-*.rpm

```

创建用户 slurm

```bash
adduser slurm

```

创建配置文件（非常关键）

```bash
mkdir -p /etc/slurm
touch /etc/slurm/slurm.conf 
```

etc中slurm.conf文件里面的配置信息来自于[https://slurm.schedmd.com/configurator.html](https://slurm.schedmd.com/configurator.html) 生成，需要配置如下选项

- SlurmctldHost: 信息来自于 `hostname -f `

- NodeName: 信息来自于 `hostname -f`, 只不过是子节点的服务器信息，如果只有单个主机，那么同上

- ComputeNodeAddress: 计算节点的IP地址，仅有单个节点时，信息为空

- PartitionName: 任务分配名，改成batch

- CPUs: 设置为空

- CoresPerSocket: 实际的物理CPU数，例如96

- ThreadsPerCore: 如果超线程，设置为2

- RealMemory: 服务器内存大小，单位为Mb

- SlurmUser: slurm要求有一个专门的用户，

- StateSaveLocation: 一定要改成 /var/spool/slurmd， 否则会出现权限问题

最后还需要增加一行 ` CgroupMountpoint=/sys/fs/cgroup `

启动 slurmctld, slurmd 的守护进程(deamon)

```bash
# 控制节点
systemctl enable slurmctld
systemctl start slurmctld
systemctl status slurmctld 
# 计算节点 
systemctl enable slurmd
systemctl start slurmd
systemctl status slurmd  
```

### 方法2: 通过OpenHPC仓库

### 测试安装

安装结果后，我们创建一个 test.sbatch, 信息如下，用于测试

```bash
#!/bin/bash
#SBATCH -J test # Job name
#SBATCH -o job.%j.out # Name of stdout output file (%j expands to %jobId)
#SBATCH -N 1 # Total number of nodes requested
#SBATCH -n 2 # Total number of mpi tasks #requested
#SBATCH -t 01:30:00 # Run time (hh:mm:ss) - 1.5 hours
# Launch MPI-based executable
echo "Test output from Slurm Testjob"
NODEFILE=`generate_pbs_nodefile`
cat $NODEFILE
sleep 20
```

递交任务

```bash
sbatch ./test.sbatch 
# Submitted batch job 2

```

查看状态

```bash
squeue
```

如果能输出一个job.X.out 文件，说明我们的SLURM已经配置成功。

## 可能报错和解决方案

使用 `rpm --install`的时候可能会遇到如下的报错。这表示你需要安装perl的Switch模块

```bahs
error: Failed dependencies:
  perl(Switch) is needed by slurm-openlava-19.05.8-1.el7.x86_64
  perl(Switch) is needed by slurm-torque-19.05.8-1.el7.x86_64

```

启动 slurmd的deamon失败

```bash
# systemctl start slurmd
Job for slurmd.service failed because the control process exited with error code. 
See "systemctl status slurmd.service" and "journalctl -xe" for details.
```

按照提示运行 `systemctl status slurmd.service` 发现error信息如下

```bash
error: Node configuration differs from hardware: Procs=1:192(hw) Boards=1:1(hw) SocketsPerBoard=1:4(hw) ...e=1:2(hw)
error: cgroup namespace 'freezer' not mounted. aborting

```

第一个error原因是在[https://slurm.schedmd.com/configurator.html](https://slurm.schedmd.com/configurator.html) 填写 "Compute Machines" 的硬件信息出现错误

第二个error原因是配置文件的默认配置表现不佳，需要做如下替换

```bash
echo CgroupMountpoint=/sys/fs/cgroup >> /etc/slurm/cgroup.conf
```

> 参考: [https://stackoverflow.com/questions/62641323/error-cgroup-namespace-freezer-not-mounted-aborting](https://stackoverflow.com/questions/62641323/error-cgroup-namespace-freezer-not-mounted-aborting)

## 参考资料

配置slurm: [https://slurm.schedmd.com/configurator.html](https://slurm.schedmd.com/configurator.html)

单节点slurm: [http://docs.nanomatch.de/technical/SimStackRequirements/SingleNodeSlurm.html](http://docs.nanomatch.de/technical/SimStackRequirements/SingleNodeSlurm.html)

munge配置:[https://github.com/dun/munge/wiki/Installation-Guide](https://github.com/dun/munge/wiki/Installation-Guide)

Slurm安装与使用: [http://wiki.casjc.com/?p=378](http://wiki.casjc.com/?p=378)

