---
title: 如何用WGDI进行共线性分析(三)
date: 2021-02-23 17:18:00.31
updated: 2021-02-24 10:49:11.06
url: /archives/comparative-genome-analysis-with-wgdi-3
categories: 基因组学
tags: 流程工具
---

在上一篇教程的最后，我写道「我们能够更加直观的观察到2个peak，基本确定2个多倍化事件」。这里就引出了一个问题，为什么我们可以根据peak来推断多倍化时间？在Lynch和Conery在2000年发表在Science的论文中，他们证明了小规模基因复制的Ks分布是L型，而在L型分布背景上叠加的峰则是来自于演化历史中某个突然的大规模复制事件。例如下图a中，实线是小规模复制的L型分布（呈指数分布, exponential distribution), 最初的峰可能是近期的串联复制引起，随着时间推移基因丢失，形成一个向下的坡。另一条虚线中的峰（呈正态分布， normal distribution）则是由全基因组复制引起。而b-h是一些物种的ks分布。

> 注，Lynch和Conery这两人是最早一批研究真核基因组的复制基因的总体保留和丢失情况的研究者。

![Kevin Vanneste et al. 2012 Molecular Biology and Evolution](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/02/image-fd0dc6527c7e416d9e5ba7b7f033a66f.png)

也就是说，我们可以通过对Ks频率分布图的观察来判断物种在历史上发生的基因组复制事件。同时又因为WGD在Ks频率分布图中表现为正态分布，那么我们可以通过对峰进行拟合来得到WGD事件所对应的Ks值。

使用WGDI对峰进行拟合需要注意的是，**它一次只能拟合一个峰**，因此我们需要先通过之前Ks可视化中所用到kspeaks模块(-kp)来过滤，接着用**PeaksFit**(-pf)模块对峰进行拟合得到模型参数，最后用KsFigures(-kf)将拟合结果绘制到一张图上。

过滤的目的是筛选出同一个WGD事件所形成的共线性区块，这可以通过设置kspeaks模块中和过滤有关的参数来实现

- pvalue: 共线性区块的显著性
- tandem: 是否过滤串联重复基因
- block_length: 共线性区块的基因对的数目
- ks_area: 根据 ks_median 筛选给定区间的共线性区块
- multiple: 选择homo的哪列作为筛选标准
- homo: 根据homo，筛选给定区间内的共线性区块对

其中比较关键的是ks_area 和homo,mutiple, 前者确定ks的大致区间，后者则是更精细地筛选共线性区块。block_length可以过滤掉一些比较小的共线性区块，毕竟基因数越多，就越不可能是随机事件所产生的共线性。

通过观察之前Ks频率分布图，我们大致确定两个峰所对应的区间分别是0~1.25，1.25~2.5。

![figure 1](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/02/image-b0c9264f436447709adeca16093f7dd9.png)

我们分别用 `wgdi -kp \? > peak1.conf`和 `wgdi -kp \? > peak2.conf`创建两个文件，主要修改其中的ks_area, homo都设置为-1,1 表示不筛选

```Ini
# peak1.conf
[kspeaks]
blockinfo  =  ath_block_information.csv
pvalue  =  0.05
tandem  =  true
block_length  =  5
ks_area  =  0,1.25
multiple  =  1
homo  =  -1,1
fontsize  =  9
area  =  0,2
figsize  =  10,6.18
savefig  =  ath.ks_peak1.distri.pdf
savefile  =  ath.ks_peak1.distri.csv

# peak2
[kspeaks]
blockinfo  =  ath_block_information.csv
pvalue  =  0.05
tandem  =  true
block_length  =  5
ks_area  =  1.25,2.5
multiple  =  1
homo  =  -1,1
fontsize  =  9
area  =  0,2
figsize  =  10,6.18
savefig  =  ath.ks_peak2.distri.pdf
savefile  =  ath.ks_peak2.distri.csv
 
```

运行 `wgdi -kp peak1.conf`  和 `wgdi -kp peak.conf`输出结果。下图中左是peak1, 右是peak2.

![figure 2 ](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/02/image-a7bf68cd528647cba17c1bcbb141ac41.png)

得到这个图之后，我们还可以进一步绘制Ks点阵图， 使用 `wgdi -bk >>peak1.conf`  和 `wgdi -bk >> peak2.conf`增加配置信息，然后对其进行修改

```Ini
# peak1
[blockks]
lens1  =  ath.len
lens2  =  ath.len
genome1_name  =  Athaliana
genome2_name  =  Athaliana
blockinfo  =  ath.ks_peak1.distri.csv
pvalue  =  0.05
tandem  =  true
tandem_length  =  200
markersize  =  1
area  =  0,2
block_length  =  5
figsize  =  8,8
savefig  =  ath.ks_peak1_dotplot.pdf
# peak2
[blockks]
lens1  =  ath.len
lens2  =  ath.len
genome1_name  =  Athaliana
genome2_name  =  Athaliana
blockinfo  =  ath.ks_peak2.distri.csv
pvalue  =  0.05
tandem  =  true
tandem_length  =  200
markersize  =  1
area  =  0,2
block_length  =  5
figsize  =  8,8
savefig  =  ath.ks_peak2_dotplot.pdf

```

运行 `wgdi -bk peak1.conf`，和 `wgdi -bk peak2.conf` 输出结果。下图中，左边是peak1，右边是peak2.

![figure 3](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/02/image-fd1f77c2359c4879869ad5865dc22c28.png)

由于我们设置的是 homo=-1,1 因此只是根据共线性块的Ks中位数筛选不同的peak，如果需要更加精细的挑选给定Ks中位数范围内的共线性区块，可以尝试修改homo值再观察结果。

> homo值表示的共线性区块中基因对的匹配情况，越接近于1，越有可能是近期WGD事件得到的基因对，越接近-1，越有可能是古老WGD事件得到的基因对。

接着，我们将用blockks步骤输出的ath.ks_peak1.distri.csv和tah.ks_peak2.distri.csv 作为peakfit模块的输入，用于拟合。

分别用 `wgdi -pf \? >> peak1.conf` 和 `wgdi -pf \? >> peak2.conf` 创建配置文件, 主要修改blockinfo

```ini
# peak1
[peaksfit]
blockinfo  =  ath.ks_peak1.distri.csv
mode  =  median
bins_number  =  200
ks_area  =  0,10
fontsize  =  9
area  =  0,3
figsize  =  10,6.18
savefig  =  ath.ks_peak1_peaksfit.pdf 

# peak2
[peaksfit]
blockinfo  =  ath.ks_peak2.distri.csv
mode  =  median
bins_number  =  200
ks_area  =  0,10
fontsize  =  9
area  =  0,3
figsize  =  10,6.18
savefig  =  ath.ks_peak2_peaksfit.pdf 
```

分别运行 `wgdi -pf peak1.conf` 和 `wgdi -pf peak2.conf`, 会得到拟合的图，以及拟合得到参数

peak1拟合参数

```Bash
R-square: 0.9804748699350867
The gaussian fitting curve parameters are :
5.02360835744403  |  0.8319599832352784  |  0.10382203381206191
```

peak2拟合参数

```Bash
R-square: 0.9188261142099129
The gaussian fitting curve parameters are :
2.084853812322066  |  1.8332872127128195  |  0.2506813629824027
```

拟合参数将会用于后续的 ksfigure 模块。

我们新建一个 all_ks.csv文件, 该文件的第一行为标题行，第二行以后为数据行。一共有4+3n列，其中第一列是样本信息，第2-3列对应线条的属性，后面都是拟合参数

```csv
 ,color,linewidth,linestyle,,,,,,
 ath_ath,green,1,--,5.02360835744403,0.8319599832352784,0.10382203381206191,2.084853812322066,1.8332872127128195,0.2506813629824027
```

注意，在终端里面编辑时，需要记住linestyle后面的逗号数量是 3n个, 其中n是最大peak数，例如拟南芥有2个peak，那么就得有6个逗号。

假如，我们之前用wgdi拟合过其他物种的peak，也得到了一些拟合参数，那么添加到这个文件中，例如

```csv
 ,color,linewidth,linestyle,,,,,,
ath_ath,red,1,--,5.02360835744403,0.8319599832352784,0.10382203381206191,2.084853812322066,1.8332872127128195,0.2506813629824027
vvi_vvi,yellow,1,-,3.00367275,1.288717936,0.177816426
vvi_oin_gamma,orange,2,-,1.910418336,1.328469514,0.262257112
vvi_oin,green,1,-,4.948194212,0.882608858,0.10426
```

最后我们用 `wgdi -kf >> ath.conf` , 创建配置文件并进行修改

```ini
[ksfigure]
ksfit = all_ks.csv
labelfontsize = 15
legendfontsize = 15
xlabel = none
ylabel = none
title = none
area = 0,4
figsize = 10,6.18
savefig =  all_ks.pdf
```

运行 `wgdi -kf  ath.conf` 得到下图

![figure 4](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/02/image-b70bff05f977435c9d22b76bd98c69ea.png)

最后，对这篇教程进行总结。

1. 我们依次使用 wgdi -kp 过滤共线性区块， wgdi -pf 对峰进行拟合， wgdi -kf 展示最后的拟合结果
2. 共线性区块的过滤是一个精细的工作，可以尝试设置不同的homo范围观察峰的变化来确定参数是否合理

