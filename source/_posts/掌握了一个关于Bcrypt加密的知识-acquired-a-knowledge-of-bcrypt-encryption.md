---
title: 掌握了一个关于Bcrypt加密的知识
date: 2023-05-31 14:12:07.55
updated: 2023-05-31 14:12:07.55
url: /archives/acquired-a-knowledge-of-bcrypt-encryption
categories: 其他
tags: 源码解读 | hash
---

> 我发现，人的表达欲是看场景的，当涉及到我最感兴趣的技术类的话题时，我其实就知道该说什么了。

最近在搞一个简单的后台管理的源码的二次开发，原项目中，不允许直接在后台添加用户，只能通过注册，或者前期初始化的时候实现。最终，在GPT的协助下，我完成了这个工作。

首先，通过阅读源码，我发现项目中的密码不是明文存放（这个年头了，不会真有人还明文吧），是经过加密的，加密代码如下

```js
const hashedPassword = bcrypt.hashSync(password, 10);
```

于是，我就让GPT根据这行命令给我写了一个脚本，结果如下

```js
// encryptPassword.js 
const bcrypt = require('bcryptjs');

function encryptPassword(password) {
  const hashedPassword = bcrypt.hashSync(password, 10);
  console.log("Hashed Password: ", hashedPassword);
  return hashedPassword;
}

// Get the password from command line arguments
const password = process.argv[2];
encryptPassword(password);
```

然后，我就发现了一个问题，那就是两次哈希的结果，居然是不一样的

```bash
$ node encryptPassword.js 123456
Hashed Password:  $2a$10$Xp6RIApHxvQJsFT2FrzNUeKfyJoHOlcSUdqHD2tU/YnnYCSuMniVK
$ node encryptPassword.js 123456
Hashed Password:  $2a$10$rDs3aQs221In0FsastR7PuarxwruyKOQ3Qd9T6p6dO81z/AUstz8q
```

我第一次遇到这个情况，于是我把结果返回给了了GPT，它给的回复如下：

 两次加密结果不一样是因为 bcryptjs 使用了盐（salt）的概念来增加密码的安全性。盐是一个随机生成的值，与密码混合后进行哈希运算，生成最终的加密结果。

每次调用 `bcrypt.hashSync` 方法时，bcryptjs 会生成一个新的随机盐，然后将该盐与密码混合后进行哈希运算。因此，即使输入的密码相同，由于使用了不同的盐，最终生成的加密结果也会不同。

这种做法是为了增加密码的安全性，即使两个用户使用相同的密码，最终的加密结果也是不同的。这样做可以避免通过比较加密结果来判断用户是否使用相同的密码。

如果需要在后续验证密码时，可以使用 `bcrypt.compareSync` 方法将输入的密码与存储的加密结果进行比较，bcryptjs 内部会自动提取盐并进行相应的比较验证。

原来如此，于是我让GPT写一个测试脚本，让我来测试下，是不是真的如此

```js
//comparePasswords.js
const bcrypt = require('bcryptjs');

function comparePasswords(password, hashedPassword) {
  const isMatch = bcrypt.compareSync(password, hashedPassword);

  if (isMatch) {
    console.log("Password matches!");
  } else {
    console.log("Password does not match!");
  }
}

// Get the user input and hashed password from command line arguments
const userInput = process.argv[2];
const hashedPassword = process.argv[3];
comparePasswords(userInput, hashedPassword);
```

我运行之后，发现两个不match啊，GPT你是不是骗我

```bash
$ node comparePasswords.js 123456 $2b$12$RJhc7WSsbRhkg2QxZuDG6.qgPXhFlrbvfOtw0SWDHkElgt1thSPD2
Password does not match
```

当我把这段内容给GPT时，它敏锐的观察到，我的输入有问题，它的回复如下

在你的示例中，你在命令行中传递了哈希值作为一个参数。这在许多命令行解释器中可能会引起问题，因为哈希字符串中的`$`符号可能会被视为变量的开头，然后尝试去解析它。结果就是，你的程序可能没有接收到完整的哈希字符串，而是只接收到了一部分。

你可以尝试用引号将哈希字符串括起来，以避免命令行解释器对`$`符号进行解析。

你的命令应该这样写：

```bash
node comparePasswords.js 123456 "$2b$12$RJhc7WSsbRhkg2QxZuDG6.qgPXhFlrbvfOtw0SWDHkElgt1thSPD2"
```

在这个版本中，哈希字符串被括在双引号中，这会防止命令行解释器将`$`符号视为变量的开头。


看来可能是因为`$`导致的。不过遗憾的是，我运行上述的代码还是失败了。因为`$`在双引号内还是会作为变量的开头，使用单引号才行，因此正确形式是

```bash
node comparePasswords.js 123456 '$2b$12$RJhc7WSsbRhkg2QxZuDG6.qgPXhFlrbvfOtw0SWDHkElgt1thSPD2'
```

不管如何，我在GPT的帮助下，了解了一个新的库，Bcrypt，还知道了它能够保证每次加密的结果都是随机的，但又保证了，同一个密码能够正确的验证。比如说我让GPT写了一个python版本，虽然输出的结果很不一样，但是依旧能够通过验证。
