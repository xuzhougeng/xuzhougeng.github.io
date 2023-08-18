---
title: 如何给bioconda贡献recipes
date: 2021-01-20 20:11:49.962
updated: 2021-01-20 20:55:38.648
url: /archives/how-to-contribute-to-bioconda
categories: Python
tags: 流程工具
---

bioconda是conda中专门用来管理生物信息相关软件的频道。

绝大部分的用户都是利用bioconda安装软件，绝大部分的教程也都是教大家如何去使用bioconda安装软件。但是bioconda的便利性并不是从天而降的，而是无数的软件开发者持续的不断的贡献和更新一些软件的recipe (菜谱) 的结果。这篇教程就是教大家如何在一个软件没有bioconda的安装方式时，自己动手写一个菜谱。

## 准备工作

首先，我们在[https://github.com/bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes) 上fork该项目到自己的GitHub项目下，然后将其clone到本地

```Bash
git clone https://github.com/xuzhougeng/bioconda-recipes.git
cd bioconda-recipes

```

如果一不小心克隆使用了官方的项目或者我的项目，也可以通过如下命令进行修改

```Bash
git remote rm origin
git remote add origin 你的github地址
```

由于bioconda是建立在conda基础上的频道，因此我们还需要预先安装好 conda，关于conda 的安装、配置和使用可以看我上传到[哔哩哔哩的视频](https://www.bilibili.com/video/BV1s4411F761)。之后安装 conda-build, 用于创建conda软件包的框架。

```Bash
# 安装conda-build, 后续用于创建模板
conda install conda-build

```

## 实际流程

向bioconda贡献一个软件的菜谱大致分为6步:

1. 创建分支

2. 编辑recipes

3. 推送修改

4. 创建pull请求

5. 删除你的分支

6. 安装你的包

## 第一步: 创建分支

我们需要先切换到主分支(master), 并保证你的主分支处于最新状态

```Bash
git checkout master
git pull upstream master
git push origin master
```

接下来，我们创建新的分支，用于处理当前的菜谱，避免影响主分支

```Bash
git checkout -b update_My_recipes
```

我建议这里的 `update_My_recipes` 修改成具体的软件名，这样更有标识性。

## 第二步: 编辑recipe

我们使用 `conda skeleton` 创建一个软件包框架， 这里pypi指的是软件托管在pypi上，wgdi则是软件名。

```Bash
# bioconda-recipes目录下
cd recipes
conda skeleton pypi wgdi 
```

因为wgdi存在一些依赖环境，如PAML, MAFFT, MUSCLE, PAL2NAL,  所以我们需要修改 `wgdi/meta.yaml` 文件，在其中加入这些软件。

```Bash
  run:
    - biopython
    - matplotlib
    - numpy
    - pandas >=1.1.0
    - python 
    - scipy
    - paml
    - mafft
    - muscle
    - pal2nal

```

事实上，conda skeleton 只是做了一小部分事情，我目前的meta.yaml 还不完整，直到后续的递交中修改了数十次才知道出问题的地方以及如何修改。

因此建议后还需要查阅 [https://docs.conda.io/projects/conda-build/en/latest/resources/define-metadata.html](https://docs.conda.io/projects/conda-build/en/latest/resources/define-metadata.html) 了解meta.yaml的更多知识。

## 第三步: 向GitHub上推送更新

这里就是用git命令将我们的修改推送到GitHub上

```Bash
# bioconda-recipes目录下
git add recipes/wgdi
git commit -m "add wgdi recipes"
git push origin update_My_recipes
```

## 第四步: 创建pull 请求

此时，让我们回到GitHub 中我们fork的bioconda-recipes页面。我们可以看到该修改已经被推送到GitHub上。


![compare & pull request](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/01/image-784347f793c94193b23914e59dc625fb.png)

我们点击 Compare & pull request, 接着就会弹出一个页面，介绍如何书写提交申请。


![Add title](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/01/image-6b3654a9dff44ff5b82cf828e96e794e.png)

这里写了4点要求，分别为

1. 在PR中加入Add或者update 作为用于说明目的，必须是第一个词

2. 如果该项目和生物学领域无关，建议换个地方。

3. 一旦通过了审核，需要在后续的issue中加入 `@BiocondaBot please add label`, 具体区别在后面的命令说明中有提到。

4. 如果有问题，你可以到Gitter或者在评论中用 `@bioconda/core` 提出自己的问题。

当你提出PR之后，bioconda 的构建系统就会测试我们提交的菜谱。如果通过了测试，我们就可以跳转到下一步，否则需要重复第二~四步，直到通过为止。

我们发现在提交之后，构建系统就提示了一个问题。


![lint problem](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/01/image-782631bcdea541fbaedca10e2ace568f.png)

点开Details，会出现如下的界面。我们可以点击View more details on CircleCI Checks 去聊了解更加详细的信息。


![view problem](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/01/image-490eb0c1569043ee89bb8b87e088b267.png)

lint问题的修改建议可在bioconda的对应页面为 [https://bioconda.github.io/contributor/linting.html](https://bioconda.github.io/contributor/linting.html) 进行查找。我把自己遇到的问题放在文末，毕竟是第一次，对于自动化测试系统并不了解，所以我反复修改了十多次。

当我们的菜谱最终通过了自动化测试后，接下来就可以在我们的issue中留言 `@bioconda-bot add label` ，让机器人增加标签，然后等待管理员review同意你的PR，当然也需要你进行修改才行。

![add label](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2021/01/image-eaa78ee47eff4337b333aa22c49e6cbd.png)

## 第五步: 删除分支

当我们的Pull Request被合并之后，我们就可以在GitHub上删除我们的分支，或者在本地删除

```Bash
# Delete local branch
git branch -D update_My_recipes
# Delete branch in your fork via the remote named "origin"
git push origin -d update_My_recipes
```

## 第六步: 安装你的包

经过前面五步之后，你需要再等待一会，才能让你的包被加入到bioconda频道中。一旦被加入到bioconda频道，那么所有人都可以非常容易的通过conda来安装你得工具了。

```Bash
conda install -c bioconda wgdi
```

## 其他

在我配置 wgdi包的时候，遇到了如下几种错误

```Bash
ERROR: recipes/wgdi/meta.yaml:0: uses_matplotlib: The recipe uses `matplotlib`, but `matplotlib-base` is recommended
ERROR: recipes/wgdi/meta.yaml:16: should_be_noarch_generic: The recipe should be build as `noarch`
ERROR: recipes/wgdi/meta.yaml:4: folder_and_package_name_must_match: The recipe folder and package name do not match.
```

根据我的错误，我在原来的yaml文件中做了如下的更改

- uses_matplotlib: 将原来的matplotlib 替换成matplobtlib-base

- should_be_noarch_generic:  在build下增加 `noarch: generic` （这一步有问题，请继续往后阅读找到正确答案，此处保留时为了记录我的思考过程）

- folder_and_package_name_must_match: 在配置文件开头增加 `{% set name = "wgdi" %}`

之后在test-linu时，我遇到了 `ModuleNotFoundError: No module named 'wgdi'`报错。这个问题和之前的 should_be_noarch_generic密切相关。noarch指的是平台无关，也就是这个软件不需要额外的编译就可以使用。一般这类软件通常是Java编译后的包，而非Python包。我最初根据提示的 `should_be_noarch_generic` 增加的 `noarch: generic` 实际并不正确，实际应该是 `noarch: python`, 也就是表示我们包依赖的Python平台。

## 参考资料:

- meta.yaml: [https://docs.conda.io/projects/conda-build/en/latest/resources/define-metadata.html](https://docs.conda.io/projects/conda-build/en/latest/resources/define-metadata.html)

- lint: [https://bioconda.github.io/contributor/linting.html](https://bioconda.github.io/contributor/linting.html) 



