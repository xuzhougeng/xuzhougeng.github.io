---
title: Singularity和MPI应用
date: 2020-06-22 01:39:27.294
updated: 2020-06-22 01:39:27.294
url: /archives/singularity-and-mpi-applications
categories: Linux
tags: 软件安装
---


MPI(Message Passin Interface)广泛应用于高性能服务器中，可用于单系统多节点或者多个计算平台间通讯，目前主流的两个开源软件分别是OpenMPI和MPICH。Singularity同时支持这两个开源工具，本篇教程介绍如何在Singularity容器开发和运行MPI程序。

在Singularity容器是执行MPI程序最常见的方法就是使用宿主机中的MPI。由于我们会同时使用宿主机的MPI和容器中的MPI，因此这被称为Host MPI或者混合模型。考虑到大部分高性能计算平台并不都支持宿主机和容器之间的文件系统共享，因此我们不讨论如何将存储盘挂载到容器，从而使用容器中宿主机MPI。

混合途径的基本想法就是，当你要执行SIngularity容器的MPI代码时，你会使用`mpiexec`类似命令去调用`singularity`命令。容器外部的MPI进程会和容器内的MPI进行写作，容器后的MPI代码会实例化任务。

Open MPI/Singularity的工作流程如下

1. 用户在shell中启动MPI(mpirun, mpiexec)
1. Open MPI接着调用进程管理守护进程(process management daemon, ORTED)
1. ORTED进程启动所需的Singularity容器
1. Singularity构建容器和命名空间环境
1. Singularity启动容器内的MPI应用
1. MPI进程启动，加载Open MPI库
1. Open MPI库通过进程管理接口(Process Management, PMI)连接回ORTED进程

此时在容器运行MPI就和直接宿主机运行MPI程序一样。优势在于它能够整合到Slurm等资源管理工具，并且和之前运行MPI应用一样。但是你需要保证容器中的MPI和宿主机的MPI版本兼容，也就是在构建MPI容器时保证和系统的MPI一致。

我们举个简单的例子。假设我们已经有了一个叫做`mpitest.c`MPI应用

```c
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char **argv) {
        int rc;
        int size;
        int myrank;

        rc = MPI_Init (&argc, &argv);
        if (rc != MPI_SUCCESS) {
                fprintf (stderr, "MPI_Init() failed");
                return EXIT_FAILURE;
        }

        rc = MPI_Comm_size (MPI_COMM_WORLD, &size);
        if (rc != MPI_SUCCESS) {
                fprintf (stderr, "MPI_Comm_size() failed");
                goto exit_with_error;
        }

        rc = MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
        if (rc != MPI_SUCCESS) {
                fprintf (stderr, "MPI_Comm_rank() failed");
                goto exit_with_error;
        }

        fprintf (stdout, "Hello, I am rank %d/%d", myrank, size);

        MPI_Finalize();

        return EXIT_SUCCESS;

 exit_with_error:
        MPI_Finalize();
        return EXIT_FAILURE;
}
```

> MPI只是一个库接口，里面的函数能够被其他编程语言所调用，所以支持Fortran, C, Python, R等。

接下来根据宿主机的MPI版本，编写一个定义文件。假如我们宿主机的MPI是MPICH，那么定义文件(命名为mpitest)如下

```bash
Bootstrap: docker
From: ubuntu:latest

%files
    mpitest.c /opt

%environment
    export MPICH_DIR=/opt/mpich-3.3
    export SINGULARITY_MPICH_DIR=$MPICH_DIR
    export SINGULARITYENV_APPEND_PATH=$MPICH_DIR/bin
    export SINGULAIRTYENV_APPEND_LD_LIBRARY_PATH=$MPICH_DIR/lib

%post
    echo "Installing required packages..."
    apt-get update && apt-get install -y wget git bash gcc gfortran g++ make

    # Information about the version of MPICH to use
    export MPICH_VERSION=3.3
    export MPICH_URL="http://www.mpich.org/static/downloads/$MPICH_VERSION/mpich-$MPICH_VERSION.tar.gz"
    export MPICH_DIR=/opt/mpich

    echo "Installing MPICH..."
    mkdir -p /tmp/mpich
    mkdir -p /opt
    # Download
    cd /tmp/mpich && wget -O mpich-$MPICH_VERSION.tar.gz $MPICH_URL && tar xzf mpich-$MPICH_VERSION.tar.gz
    # Compile and install
    cd /tmp/mpich/mpich-$MPICH_VERSION && ./configure --prefix=$MPICH_DIR && make install
    # Set env variables so we can compile our application
    export PATH=$MPICH_DIR/bin:$PATH
    export LD_LIBRARY_PATH=$MPICH_DIR/lib:$LD_LIBRARY_PATH
    export MANPATH=$MPICH_DIR/share/man:$MANPATH

    echo "Compiling the MPI application..."
    cd /opt && mpicc -o mpitest mpitest.c
```

如果是Open MPI,

```bash
Bootstrap: docker
From: ubuntu:latest

%files
    mpitest.c /opt

%environment
    export OMPI_DIR=/opt/ompi
    export SINGULARITY_OMPI_DIR=$OMPI_DIR
    export SINGULARITYENV_APPEND_PATH=$OMPI_DIR/bin
    export SINGULAIRTYENV_APPEND_LD_LIBRARY_PATH=$OMPI_DIR/lib

%post
    echo "Installing required packages..."
    apt-get update && apt-get install -y wget git bash gcc gfortran g++ make file

    echo "Installing Open MPI"
    export OMPI_DIR=/opt/ompi
    export OMPI_VERSION=4.0.1
    export OMPI_URL="https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-$OMPI_VERSION.tar.bz2"
    mkdir -p /tmp/ompi
    mkdir -p /opt
    # Download
    cd /tmp/ompi && wget -O openmpi-$OMPI_VERSION.tar.bz2 $OMPI_URL && tar -xjf openmpi-$OMPI_VERSION.tar.bz2
    # Compile and install
    cd /tmp/ompi/openmpi-$OMPI_VERSION && ./configure --prefix=$OMPI_DIR && make install
    # Set env variables so we can compile our application
    export PATH=$OMPI_DIR/bin:$PATH
    export LD_LIBRARY_PATH=$OMPI_DIR/lib:$LD_LIBRARY_PATH
    export MANPATH=$OMPI_DIR/share/man:$MANPATH

    echo "Compiling the MPI application..."
    cd /opt && mpicc -o mpitest mpitest.c
```

接着我们从定义文件中构建出sig文件, `singularity build mpitest.sig mpitest`

最后，执行Singularity容器的MPI应用的标准方法就是在host中运行`mpirun`命令，它会启动Singularity容器，并最终运行其中的MPI程序。

```bash
mpirun -n 4 singularity exec mpitest.sig /opt/mpitest
# 输出信息如下
# Hello, I am rank 0/4Hello, I am rank 1/4Hello, I am rank 2/4Hello, I am rank 3/4
```

如果是在SLURM这类任务提交系统中，我们需要批处理脚本来执行MPI引用，举个例子

```bash
#!/bin/bash
#SBATCH --job-name singularity-mpi
#SBATCH -N $NNODES # total number of nodes
#SBATCH --time=00:05:00 # Max execution time

mpirun -n $NP singularity exec /var/nfsshare/gvallee/mpich.sif /opt/mpitest
```

执行方式如下

```bash
sbatch my_job.sh
```

参考资料: https://sylabs.io/guides/3.3/user-guide/mpi.html