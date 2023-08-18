---
title: 基于SpringBoot和React的一个文件上传案例
date: 2023-02-22 11:12:48.103
updated: 2023-02-22 11:12:48.103
url: /archives/springboot-and-react-file-upload
categories: Java
tags: Web开发
---

利用后端的SpringBoot框架和前端的React框架实现一个文件上传的小项目，对目前的学习做一个总结。

# 环境搭建

我们先分别搭建前后端的项目环境，然后进行项目开发。

## 前端环境搭建

> 前端开发需要一些基本的HTML, CSS和Javscript的背景知识。

我们直接使用成前端比较热门的React框架进行网页搭建，还有一个热门框架是Vue。

先安装NODE环境，在[https://nodejs.org/en/](https://nodejs.org/en/)下载长期支持版，安装之后，在命令行输入node -v 确定版本。

编程环境选择VSCode，同时安装了几个插件：

- IntelliJ IDEA Keybindings
- Prettier
- Prettier ESlint

之后，我们使用[umi](https://umijs.org)作为项目的脚手架，搭建项目

```Bash
mkdir file_upload_frontend && cd file_upload_frontend
pnpm dlx create-umi@latest
# 选择simple App + pnpm + taobao
# 安装ant-design
pnpm i antd @ant-design/icons

```

安装结束后，启动项目，访问Network里的地址，就确定项目是跑上了。

```Bash
pnpm dev
....
        ╔════════════════════════════════════════════════════╗
        ║ App listening at:                                  ║
        ║  >   Local: http://localhost:8000                  ║
ready - ║  > Network: http://192.168.50.97:8000              ║
        ║                                                    ║
        ║ Now you can open browser with the above addresses↑ ║
        ╚════════════════════════════════════════════════════╝
...

```

网页显示如下

![启动页面](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/02/image-778d8da29fb04b72add7704afe96479f.png)

## 后端环境搭建

后端使用Java SpringBoot实现RESTful的接口开发。(这里的背景知识要求比较多，Java的面向对象编程，反射，注解，Spring框架，SpringMVC框架）

首先，到[https://www.jetbrains.com/idea/download/](https://www.jetbrains.com/idea/download/)里下载一个社区版的IntelliJ IDEA。相比于商业版，社区办缺少了Spring初始化、SQL数据库连接、HTTP客户端等功能，但是也是能用。

接着，在[https://start.spring.io](https://start.spring.io)中创建一个SpringBoot项目，配置如下

![image.png](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/02/image-4f5fa9ed7847457899a4eaf59b2c2932.png)

得到一个fileupload.zip文件，解压缩后就得到了我们项目的基本结构

```text
fileupload
├── HELP.md
├── mvnw
├── mvnw.cmd
├── pom.xml
└── src
    ├── main
    │   ├── java
    │   │   └── top
    │   │       └── xuzhougeng
    │   │           └── fileupload
    │   │               └── FileuploadApplication.java
    │   └── resources
    │       ├── application.properties
    │       ├── static
    │       └── templates
    └── test
        └── java
            └── top
                └── xuzhougeng
                    └── fileupload
                        └── FileuploadApplicationTests.java
```

利用IDEA打开这个项目，把resources下的application.properties改成application.yaml，编辑内容如下

```YAML
server:
  port: 8070
```

然后点击IDEA的运行键，如果能正常访问localhost:8070就意味着后端环境配置成功

![springboot](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/02/image-b39e7be42bd542b19113f956bd19ca37.png)

# 代码编写

## 前端基础代码

前端的工作非常简单，只需要去Ant Design的官方找一个合适的模版。

打开[https://ant-design.gitee.io/components/overview-cn/](https://ant-design.gitee.io/components/overview-cn/)，找到Upload上传组件块，选择一个顺眼的，复制它的代码

![上传模块](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/02/image-0057b758c1124b98a8c9af9413f006a0.png)

之后在项目中建立一个src/pages/upload.tsx文件，把代码粘贴进去

修改 .umirc.tsx 中的路由配置

```React TSX
import { defineConfig } from "umi";

export default defineConfig({
  routes: [
    { path: "/", component: "index" },
    { path: "/docs", component: "docs" },
    { path: "/upload", component: "upload"}
  ],
  npmClient: 'pnpm',
});

```

打开[http://localhost:8000/upload](http://localhost:8000/upload) 如果看到里面内容跟复制的内容一样，就算成功。

接下来，编辑upload.tsx文件，做一些基本的修改

```React TSX
  multiple: false, //只允许传一个文件
  action: 'http://localhost:8070/files', //设置上传的站点
```

multiple设置为false，用于限制单次上传的文件数，action设置上传的目标地址。

这样子，前端就算完工了。

## 后端代码

后端文件上传代码参考[https://www.bezkoder.com/spring-boot-file-upload/](https://www.bezkoder.com/spring-boot-file-upload/)

首先，我们写服务层的代码，实现两个功能

1. 初始化文件上传后存放的目录
2. 保存文件

在top.xuzhougeng.fileupload建立一个软件包, service。

然后在service中创建一个接口FileStorageService，里面就定义两个接口函数

```Java
// top/xuzhougeng/fileupload/service/FileStorageService.java
package top.xuzhougeng.fileupload.service;
import org.springframework.web.multipart.MultipartFile;

public interface FileStorageService {

    public void init();
    public void save(MultipartFile file);
}

```

通过在FileStorageService上用快捷键option + enter，实现接口。

![生成方法](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/02/image-060d043e340849a192aea10ad2193439.png)

实现类存放在service的impl目录下

![实现类](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2023/02/image-f7afc15de3bf4445826a052c19f8eb2f.png)

代码如下

```Java
// top/xuzhougeng/fileupload/service/FileStorageServiceImpl.java
package top.xuzhougeng.fileupload.service.impl;

import org.springframework.stereotype.Service;
import org.springframework.web.multipart.MultipartFile;
import top.xuzhougeng.fileupload.service.FileStorageService;

import java.io.IOException;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;


@Service
public class FileStorageServiceImpl implements FileStorageService {

    private final Path root = Paths.get("uploads");
    @Override
    public void init() {
        try {
            Files.createDirectories(root);
        } catch (IOException e){
            throw new RuntimeException("无法初始化upload的上传目录");
        }

    }

    @Override
    public void save(MultipartFile file) {
        try {
            Files.copy(file.getInputStream(),
                    this.root.resolve(file.getOriginalFilename()));
        } catch (Exception e){
            if (e instanceof FileAlreadyExistsException){
                throw new RuntimeException("文件已存在");
            }

            throw new RuntimeException(e.getMessage());
        }
    }
}


```

初始化init调用Files模块创建文件夹。save函数的参数类型为MultipartFile，是SpringMVC用来简化文件上传操作的工作类。通过调用Files模块，将上传的文件流复制到指定的目录下。

服务层搞定之后，就是控制层，创建一个软件包controller，并在目录下建立一个FileController类。

```Java
package top.xuzhougeng.fileupload.controller;


import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;
import top.xuzhougeng.fileupload.service.FileStorageService;


@RestController
@CrossOrigin(value = {"http://192.168.50.97:8000"})
public class FileController {

    private FileStorageService fileStorageService;
    @Autowired
    public void setFileStorageService(FileStorageService fileStorageService) {
        this.fileStorageService = fileStorageService;
    }

    @PostMapping("/files")
    public String fileSave(@RequestParam("file") MultipartFile file){
        System.out.println("loading file");
        try {
            fileStorageService.save(file);
            return "Upload successfully";
        } catch (Exception e){
            return "Upload failed";
        }

    }
}

```

它的功能非常简单，就是相应POST请求，然后调用服务层的fileStorageService的save方法进行文件保存。

最后还需要修改启动函数，使其实现CommandLineRunner，覆写他的 run方法，使其在运行前调用FileStorageService的init方法。

```Java
@SpringBootApplication
public class FileuploadApplication implements CommandLineRunner {

  @Resource
  private FileStorageService fileStorageService;

  @Override
  public void run(String... args) throws Exception {
    fileStorageService.init();
  }

  public static void main(String[] args) {
    SpringApplication.run(FileuploadApplication.class, args);
  }
  

}
```

最后启动后端，测试下文件上传。
