---
title: 无root权限解决编译时的依赖问题
date: 2019-09-23 11:02:01.221
updated: 2019-09-23 11:02:01.221
url: /archives/Compile-Software-in-Linux-Without-Root
categories: Linux
tags: 软件安装
---

如果你拥有最高权限，如果你只管理一台服务器，那么系统自带的包管理工具就帮你解决了所有问题。但是真实世界没有那么美好，所以我花了那么久时间去学习如何从源码开始编译一个软件。

**环境**为CentOS Linux release 7.4.1708 (Core), Linux内核version 3.10.0-693.el7.x86\_64， GCC版本为4.8.5 20150623 (Red Hat 4.8.5-16) (GCC)，

## Linux的编译体系

无管理员权限编译的常规三部曲是`./configure --prefix=$HOME/usr && make && make install`，其中最重要的一步就是`configure`，它所做的任务如下

- 检查GCC版本以及是否安装了编译所需工具
- 如果需要头文件，则默认去`/usr/include`查找
- 如果涉及到动态编译库，则默认去`/usr/lib`和`/usr/lib64`查找. 注：`lib`的函数库仅用于开机时用,提供给/bin和/sbin.

那为何需要配置？配置主要解决软件开发和软件实际安装时平台不同所导致的问题。

首先，一个C/C++工程不可能只用到标准库，很多已有的轮子就不需要重复制造。其次，由于很多软件都重复用到相同的依赖库，那么如果把这些依赖库独立成单独的模块，在调用的时候加载，也能节省空间。早期为了适配多个平台，开发人员需要手写条件语句来检查环境依赖，后来GNU专门开发了Autotools辅助构建源码编译所需要的关键文件。

![Autotools](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/9/2013053-36054cdebaf2ec9d-d308d285e8734e5eb72f778be7e33c98.png)

## 编译环境变量

用`./configure -h`查看帮助文档的时候，会提示几个和编译相关非常重要的环境变量。

```shell
# 编译器
CC          指定C编译器(compiler command)路径
CXX         指定C++编译器
# 编译器选项
CFLAGS      用于C编译器的选项
CXXFLAGS    用于C++编译器的选项
LDFLAGS     链接相关选项,如果你有自定义的函数库(lib dir)，即可以用 -L<lib dir>指定
# 预编译器
CXXCPP      C++ 预处理器(preprocessor)
CPP         C 预处理器(preprocessor)
# 预编译器选项
CPPFLAGS    C/C++预处理器选项, 如果你自定义的头文件，可以用-I<include dir>
```

Makfile规则中的编译命令通常遵守如下规范：

1,首先从源代码生成目标文件( **预处理,编译,汇编** )，"-c"选项表示不执行链接步骤；

```shell
$(CC) $(CPPFLAGS) $(CFLAGS) example.c   -c   -o example.o
```

2,然后将目标文件链接为最终的结果( **链接** )，"-o"选项用于指定输出文件的名字。

```shell
$(CC) $(LDFLAGS) example.o   -o example
```

这些只是约定俗成的惯例，所以有些人会“随性而为”，你也拿他没有办法。尽管将源代码编译为二进制文件的四个步骤由不同的程序(cpp,gcc/g++,as,ld)完成，但是事实上 cpp, as, ld 都是由 gcc/g++ 进行间接调用的。换句话说，**控制了 gcc/g++ 就等于控制了所有四个步骤**。从 Makefile 规则中的编译命令可以看出，编译工具的行为全靠 **CC/CXX CPPFLAGS CFLAGS/CXXFLAGS LDFLAGS** 这几个变量在控制。所以控制这些变量最简单的做法是首先设置与这些 Makefile 变量同名的环境变量并将它们 export 为 **环境变量**（全局），然后运行 configure 脚本，大多数 configure 脚本会使用这同名的环境变量代替 Makefile 中的值

- CC/CXX: 指定C/C++编译所在路径，即可以选择不同的版本的编译器进行编译。
- CPPFLAGS: 这是用于预处理阶段的选项。用于添加不在标准路径`/usr/include`下的头文件。如`CPPFLAGS="-I$HOME/usr/include -I$HOME/usr/include/ncurses"`
- CFLAGS/CXXFLAGS： 指定头文件（.h文件）的路径，如：`CFLAGS=-I/usr/include -I/path/include`。同样地，安装一个包时会在安装路径下建立一个include目录，当安装过程中出现问题时，试着把以前安装的包的include目录加入到该变量中来。

> CPPFLAG和CFLAGS/CXXFLAGS这两个变量可以认为是等价关系，都用在预处理阶段，也就是都能用于指定头文件所在位置。

- LDFLAGS：gcc 等编译器会用到的一些优化参数，也可以在里面指定库文件的位置。用法：`LDFLAGS=-L/usr/lib -L/path/to/your/lib`。每安装一个包都几乎一定的会在安装目录里建立一个lib目录。如果明明安装了某个包，而安装另一个包时，它愣是说找不到，可以抒那个包的lib路径加入的LDFALGS中试一下。

有时候LDFLAGS指定-L虽然能让链接器找到库进行链接，但是运行时链接器却找不到这个库，如果要让软件运行时库文件的路径也得到扩展，那么我们需要增加这两个库给"-Wl,R"：

```shell
LDFLAGS = -L/var/xxx/lib -L/opt/mysql/lib -Wl,R/var/xxx/lib -Wl,R/opt/mysql/lib
```

如在执行./configure以前设置环境变量 `export LDFLAGS="-L/var/xxx/lib -L/opt/mysql/lib -Wl,R/var/xxx/lib -Wl,R/opt/mysql/lib"`，注意设置环境变量等号两边不可以有空格，而且要加上引号（shell的用法）。那么执行configure以后，Makefile将会设置这个选项，链接时会有这个参数，编译出来的可执行程序的库文件搜索路径就得到扩展了

除了通过以上几种环境变量为编译器提供头文件和静态和动态库文件的位置信息外，还有一种变量叫做 **PKG\_CONFIG\_PATH** , 它从`xx.pc`文件获取读取相应的环境环境。

**注意**:Linux下编译共享库时，必须加上-fPIC参数，即`export CFLAGS=" -fPIC" CXXFLAGS=" -fPIC"`否则在链接时会有错误提示.这是在编译zsh时候发现明明装了ncurse却还是不能用的共享库的坑。

> fPIC的目的是什么？共享对象可能会被不同的进程加载到不同的位置上，如果共享对象中的指令使用了绝对地址、外部模块地址，那么在共享对象被加载时就必须根据相关模块的加载位置对这个地址做调整，也就是修改这些地址，让它在对应进程中能正确访问，而被修改到的段就不能实现多进程共享一份物理内存，它们在每个进程中都必须有一份物理内存的拷贝。fPIC指令就是为了让使用到同一个共享对象的多个进程能尽可能多的共享物理内存，它背后把那些涉及到绝对地址、外部模块地址访问的地方都抽离出来，保证代码段的内容可以多进程相同，实现共享。

**动态库路径问题**: 由前面可以知道许多大型软件为了减少体积不会完全编译所有功能，而是借助于动态连接技术在运行时加载封装好的动态连接库内的功能。这就涉及一个非常重要的问题，软件如何知道动态链接库所处的位置。动态库搜索路径分两种情况，一种是编译生成可执行文件时，另外一种是运行可执行文件时。

编译生成可执行文件时，动态库的搜索路径顺序如下：

- 首先gcc会找-L选项；
- 然后再找gcc的环境变量`LIBRARY_PATH`，可以在.profile设置这个环境变量，并且可以通过选项-v查看gcc最终编译时LIBRARY_PATH的值；
- 再找内定目录: `/lib：/usr/lib：/usr/local/lib`，这些都是当初compile gcc时写在程序内的。

注意上面索顺序是**不会递归**在目录下搜索的。

生成可执行文件后，运行文件时，动态库的搜索路径顺序如下：

- 首先编译目标代码时指定的动态库搜索路径，就是用选项 `-Wl,rpath` 指定程序在运行时动态库的搜索路径，比如gcc -Wl,-rpath,include -L. -ldltest hello.c，在执行文件时会搜索路径./include；
- 环境变量`LD_LIBRARY_PATH`指定的动态库搜索路径；
- 配置文件`/etc/ld.so.conf`中指定的动态库搜索路径，即在配置文件中添加动态库的绝对路径，然后运行指令ldconfig是配置文件生效；
- 默认的动态库搜索路径`/lib:/usr/lib`。

同样上面索顺序是不会递归在目录下搜索的。通常使用动态库简单做法是：把生成的so文件拷贝到/usr/lib中，这样不管是生成可以执行文件时，还是执行程序时，都能找到需要的so文件。但是普通用户没有/usr/lib的写入权限，所有要指定`LD_LIBRARY_PATH`.ls

参考资料:

- [CFLAGS详解](http://blog.csdn.net/xinyuan510214/article/details/50457433)
- [Makefile编译选项CC与CXX/CPPFLAGS、CFLAGS与CXXFLAGS/LDFLAGS](http://blog.csdn.net/hjwang1/article/details/44497489)
- [使用gcc时头文件路径和动态链接库路径](http://blog.csdn.net/MaximusZhou/article/details/38559963)

## GCC安装(非必要)

首先让我们利用系统原来老旧的GCC编译器编译出最新版本的gcc吧，毕竟安装软件的时候，GCC的版本一定要过最低要求。

**第一步**： 下载gcc源码

```shell
mkdir -p ~/src && cd ~/src
wget https://mirrors.tuna.tsinghua.edu.cn/gnu/gcc/gcc-7.2.0/gcc-7.2.0.tar.gz
tar -zxvf gcc-7.2.0.tar.gz && cd gcc-7.2.0
ls
```

![check](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/9/2013053-c5d4c17d7e569961-e411ca5a623a4f4ca75bd3596b0b5530.png)

**第二步**， 检查系统是否已经具备前置软件, 主要是GMP，MPFR, MPC。这些软件可以到<ftp://gcc.gnu.org/pub/gcc/infrastructure/>找到，然后下载后解压缩，并移动到gcc源码文件夹下。 可以在配置的时候用`--with-gmp, --with-mpfr --with-mpc`指定具体所在路径。

```shell
cd src
# GNU Multiple precision Library
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/gmp-6.1.0.tar.bz2 \
&& tar -jxvf gmp-6.1.0.tar.bz2 && mv gmp-6.1.0 gcc-7.2.0/gmp
# isl library
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/isl-0.18.tar.bz2 \
&& tar -jxvf isl-0.18.tar.bz2 && mv isl-0.18 gcc-7.2.0/isl
# MPFR Library
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/mpfr-3.1.4.tar.bz2 \
&& tar -jxvf mpfr-3.1.4.tar.bz2 && mv mpfr-3.1.4 gcc-7.2.0/mpfr
# MPC Library
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/mpc-1.0.3.tar.gz \
&& tar -zxvf mpc-1.0.3.tar.gz && mv mpc-1.0.3 gcc-7.2.0/mpc
```

不过更加人性化的方法是在GCC源码根目录下运行`./contrib/download_prerequisites`，可以自动搞定。

**第三步**：使用`./configure`进行配置。官方**强烈**建议, 配置所在文件夹一定要和源码所在文件夹区分开，此外configure还可以配置很多参数，我的代码如下：

```shell
mkdir build && cd build
../configure\
	--prefix=$HOME/usr \ # 指定安装路径
	--disable-multilib \ # 取消32位库编译
	--enable-threads=posix \ # 使用POSIX/Unix98作为线程支持库
```

基本上这一步不会出现太多的报错，都能够顺利生成Makefile.

**第四步**： 编译. 这步有一个小技巧就是利用多核处理器进行加速，例如`make -j2` 就是双核。

这一部分很慢很慢，因为默认设置下是3个阶段的引导(3-stage bootstrap), 以保证能够编译出完整的GCC系统并且还不会出错，你可以在配置的时候用`--disable-bootstrap`进行关闭。

**第五步**： 安装。如果你编译都成功了，那么安装也不会存在问题了， `make install`.

那么我们编译的GCC和系统自带的有什么**区别**吗？

```shell
# 从头编译
$ $HOME/usr/bin/gcc -v
Using built-in specs.
COLLECT_GCC=/home/zgxu/usr/bin/gcc
COLLECT_LTO_WRAPPER=/home/zgxu/usr/libexec/gcc/x86_64-pc-linux-gnu/7.2.0/lto-wrapper
Target: x86_64-pc-linux-gnu
Configured with: ../configure --prefix=/home/zgxu/usr --disable-multilib --enable-threads=posix
Thread model: posix
gcc version 7.2.0 (GCC)
# 系统自带
$ gcc -v
Using built-in specs.
COLLECT_GCC=gcc
COLLECT_LTO_WRAPPER=/usr/libexec/gcc/x86_64-redhat-linux/4.8.5/lto-wrapper
Target: x86_64-redhat-linux
Configured with: ../configure --prefix=/usr --mandir=/usr/share/man --infodir=/usr/share/info --with-bugurl=http://bugzilla.redhat.com/bugzilla --enable-bootstrap --enable-shared --enable-threads=posix --enable-checking=release --with-system-zlib --enable-__cxa_atexit --disable-libunwind-exceptions --enable-gnu-unique-object --enable-linker-build-id --with-linker-hash-style=gnu --enable-languages=c,c++,objc,obj-c++,java,fortran,ada,go,lto --enable-plugin --enable-initfini-array --disable-libgcj --with-isl=/builddir/build/BUILD/gcc-4.8.5-20150702/obj-x86_64-redhat-linux/isl-install --with-cloog=/builddir/build/BUILD/gcc-4.8.5-20150702/obj-x86_64-redhat-linux/cloog-install --enable-gnu-indirect-function --with-tune=generic --with-arch_32=x86-64 --build=x86_64-redhat-linux
Thread model: posix
gcc version 4.8.5 20150623 (Red Hat 4.8.5-16) (GCC)
```

不谈安装路径和版本，基本上 **差别** 就是在配置这一步，而这些参数就需要仔细研究了。

一个 **错误** : 'Link tests are not allowed after GCC\_NO\_EXECUTABLES.' 后来发现是第三步没有在独立的文件下构建Makefile.

参考资料：

- installing GCC: <https://gcc.gnu.org/install/>
- linux下编译gcc6.2.0: <https://www.cnblogs.com/oloroso/p/5984985.html>

## CMake: 平台无关的编译软件

不同平台有着不同的Make工具用于编译,例如 GNU Make ，QT 的 qmake ，微软的 MS nmake，BSD Make（pmake），Makepp，等等。这些 Make 工具遵循着不同的规范和标准，所执行的 Makefile 格式也千差万别。这样就带来了一个严峻的问题：如果软件想跨平台，必须要保证能够在不同平台编译。而如果使用上面的 Make 工具，就得为每一种标准写一次 Makefile ，这将是一件让人抓狂的工作。

CMake就是针对上面问题所设计的工具：它首先允许开发者编写一种平台无关的 CMakeList.txt 文件来定制整个编译流程，然后再根据目标用户的平台进一步生成所需的本地化 Makefile 和工程文件，如 Unix 的 Makefile 或 Windows 的 Visual Studio 工程。从而做到“Write once, run everywhere”。显然，CMake 是一个比上述几种 make 更高级的编译配置工具。一些使用 CMake 作为项目架构系统的知名开源项目有 VTK、ITK、KDE、OpenCV、OSG 等.

```bash
wget https://cmake.org/files/v3.10/cmake-3.10.2.tar.gz
tar xf cmake-3.10.2.tar.gz
cd cmake-3.10.2
```

参考资料：

- <http://www.hahack.com/codes/cmake/>
- <https://www.cnblogs.com/d-blog/p/4617208.html>

## 几个必须要装的函数库

在安装之前需要先声明几个环境变量，可以直接添加在配置文件中。这都是后面遇到共享库的问题得到的经验教训。

```shell
export CFLAGS=" -fPIC"
export CXXFLAGS=" -fPIC"
export CPPFLAGS="-I$HOME/usr/include -I$HOME/usr/include/ncurses -I$HOME/usr/include/X11"
export LDFLAGS="-L$HOME/usr/lib -L$HOME/usr/lib64"
export LD_LIBRARY_PATH=$HOME/usr/lib:$HOME/usr/lib64
export PKG_CONFIG_PATH=$HOME/usr/lib/pkgconfig:$HOME/usr/share/pkgconfig
```

**ncurses**提供了一系列的函数以便使用者调用它们去生成基于文本的用户界面，许多大名鼎鼎的软件都用到了ncurses，例如vim, screen,tmux,zsh等。并且**samtools**如果需要tview可视化BAM文件，也需要这个库做支持。

```shell
wget ftp://ftp.invisible-island.net/ncurses/ncurses.tar.gz && tar -zxvf ncurses.tar.gz
./configure --enable-shared --prefix=$HOME/usr
make && make install
```

**Libevent**是一个用C语言编写的、轻量级的开源高性能事件通知库, 后续安装tmux时候需要这个依赖库。

```shell
# libevent
cd src
wget https://github.com/libevent/libevent/releases/download/release-2.1.8-stable/libevent-2.1.8-stable.tar.gz
tar -zxvf libevent-2.1.8-stable.tar.gz && cd  libevent-2.1.8
./configure prefix=$HOME/usr && make && make install
```

**bzip2, xz, zlib**: 文件压缩相关函数库，后续samtools编译时需要。

```shell
wget http://www.zlib.net/zlib-1.2.11.tar.gz
tar -zxvf zlib-1.2.11.tar.gz && cd zlib-1.2.11 && ./configure --prefix=$HOME/usr && make && make install
wget http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz
tar -zxvf bzip2-1.0.6.tar.gz && cd bzip2-1.0.6 && ./configure --prefix=$HOME/usr && make && make install
wget https://tukaani.org/xz/xz-5.2.3.tar.gz
tar -zxvf xz-5.2.3.tar.gz && cd xz-5.2.3 && ./configure --prefix=$HOME/usr && make && make install
```

**openssl, libssh2, libcurl**: 计算机之间文件传输访问相关库。其中OpenSSL是一个安全套接字层密码库，囊括主要的密码算法、常用的密钥和证书封装管理功能及SSL协议，并提供丰富的应用程序供测试或其它目的使用。libssh2是一个C 函数库，用来实现SSH2协议。libcurl主要功能就是用不同的协议连接和沟通不同的服务器.

```shell
# 安装有先后
# openssl
wget https://www.openssl.org/source/openssl-1.0.2m.tar.gz
tar -zxvf openssl-1.0.2m.tar.gz && cd openssl-1.0.2m
# 这里非常神奇的居然是config，添加shared生成动态库
./config prefix=$HOME/usr shared
make && make install
# 卸载使用 make clean
# libssh2
wget https://www.libssh2.org/download/libssh2-1.8.0.tar.gz
tar -zxvf libssh2-1.8.0.tar.gz && cd libssh2-1.8.0
./configure --with-libssl-prefix=$HOME/usr/ssl --prefix=$HOME/usr
# libcurl
wget https://curl.haxx.se/download/curl-7.56.1.tar.gz
tar -zxvf curl-7.56.1.tar.gz && cd curl-7.56.1
./configure --prefix=$HOME/usr --enable-http --enable-ftp --enable-file --enable-proxy --enable-telnet --enable-libcurl-option --enable-ipv6 --with-lib --with-ssl
```

**readline**: GNU提供用于这些命令补全、搜索历史命令、行编辑快捷键等等这些人性化的交互方式的函数库，缺少这个标准库，编译的R就缺少自动补全的功能。

```shell
wget http://ftp.gnu.org/gnu/readline/readline-7.0.tar.gz
tar -zxvf readline-7.0.tar.gz && cd readline-7.0
./configure --prefix=$HOME/usr && make && make install
```

**PCRE**: 提供和Perl5相同语法和语义正则表达式的函数库，后续安装R用到。

```shell
wget ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.41.tar.gz
tar -zxvf pcre-8.41.tar.gz && cd pcre-8.41
./configure --enable-utf --enable-pcregrep-libz --enable-pcregrep-libbz2 --prefix=$HOME/usr
```

**x11**：X11也叫做X Window系统，X Window系统 (X11或X)是一种位图显示的视窗系统,是在 Unix 和 类Unix 操作系统，以及OpenVMS上建立图形用户界面的标准工具包和协议，并可用于几乎所有已有的现代操作系统。主要是R编译的时候要用，具体用途有待开发。

x11安装比较复杂，有很多的依赖库，因此需要先安装xtrans, xextproto, xcb(lib,xcb-proto, libpthread-subs), kbproto,xproto,inputproto。但是编译很容易，仅提供下载链接

```shell
https://www.x.org/releases/X11R7.7/src/lib/xtrans-1.2.7.tar.gz
https://www.x.org/releases/X11R7.7/src/proto/xextproto-7.2.1.tar.gz
https://www.x.org/releases/X11R7.7/src/proto/kbproto-1.0.6.tar.gz
https://www.x.org/releases/X11R7.7/src/proto/xproto-7.0.23.tar.gz
https://www.x.org/releases/X11R7.7/src/proto/inputproto-2.2.tar.gz
https://www.x.org/releases/X11R7.7/src/xcb/libpthread-stubs-0.3.tar.gz
https://www.x.org/releases/X11R7.7/src/xcb/xcb-proto-1.7.1.tar.gz
https://www.x.org/releases/X11R7.7/src/xcb/libxcb-1.8.1.tar.gz
https://www.x.org/releases/X11R7.7/src/lib/libSM-1.2.1.tar.gz
https://www.x.org/releases/X11R7.7/src/lib/libICE-1.0.8.tar.gz
https://www.x.org/releases/X11R7.7/src/lib/libXt-1.1.3.tar.gz
```

相当于人工检查依赖环境，仅仅是繁琐而已，并不困难

```shell
# 安装X11
wget -4 https://www.x.org/releases/X11R7.7/src/lib/libX11-1.5.0.tar.gz
tar -zxvf libX11-1.5.0.tar.gz && cd libX11-1.5.0
./configure --prefix=$HOME/usr && make && make install
```

## 编译案例

### 安装zsh

zsh或许可以认为是最好的shell，用过zsh的人都不会想着bash了。不过zsh自定义配置，可供选择的插件以及主题实在是太多，因此一定要搭配oh-my-zsh。zsh依赖ncurses.

```shell
wget -O zsh.tar.gz https://sourceforge.net/projects/zsh/files/latest/download
tar -zxvf zsh.tar.gz && cd zsh
export CPPFLAGS="-I$HOME/usr/include/" LDFLAGS="-L$HOME/usr/lib"
../configure --prefix=$HOME/usr --enable-shared
make && make install
```

由于没有root权限，无法使用`chsh`，只能通过在`~/.bashrc`添加`exec $HOME/usr/bin/zsh -l`保证登陆的时候自动切换成zsh。其次, zsh搭配oh-my-zsh才完整, 只不过这里只能手动安装了。

```shell
# 从github上克隆oh-my-zsh
git clone git://github.com/robbyrussell/oh-my-zsh.git ~/.oh-my-zsh
# 用oh-my-zsh的zsh配置文件替代
cp ~/.oh-my-zsh/templates/zshrc.zsh-template ~/.zhsrc
# 安装一些字体, 不然一些主题会显示异常
cd src
git clone https://github.com/powerline/fonts.git --depth=1
cd fonts && ./install.sh
```

重启一下终端，后面根据需要调整配置文件里的参数。

### 编译tmux

tmux和screen类似，也是文本终端神器, 依赖于libevent和ncurses.

```shell
export CPPFLAGS="-I$HOME/usr/include -I$HOME/usr/include/ncurses"
export LDFLAGS="-L$HOME/usr/lib -L$HOME/usr/lib64"
mkdir -p src && cd src
git clone https://github.com/tmux/tmux.git
cd tmux
sh autogen.sh
./configure --prefix=$HOME/usr
make && make install
```

### 编译R语言

由于我自己编译完全版的GCC套餐，很多之前的gfortran不存在的问题也就不存在了（管理员安装了Java）。此外，R还需要gnu readline, pcre > 8.2, x11。当然这些函数包都在之前安装好了。

一些依赖库

```shell
# 安装tidyverse发现xm12需要libiconv的libiconv.so.2
https://ftp.gnu.org/pub/gnu/libiconv/libiconv-1.15.tar.gz
```

正式安装

```shell
wget https://cran.r-project.org/src/base/R-3/R-3.4.2.tar.gz
tar -zxvf R-3.4.2.tar.gz  && cd R-3.4.2/
./configure --prefix=$HOME/R
make && make install
```

![R configure](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2019/9/2013053-2bae1c5cdc39de29-b1856615b47a43b692e8da39655adcdc.png)

到此，我可以说Linux平台下即便我没有root权限，也没有多少软件包是我所不能手工编译。

