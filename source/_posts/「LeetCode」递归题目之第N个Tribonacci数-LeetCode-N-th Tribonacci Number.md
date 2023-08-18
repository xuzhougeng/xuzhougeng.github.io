---
title: 「LeetCode」递归题目之第N个Tribonacci数
date: 2019-10-22 21:21:35.374
updated: 2019-10-22 21:22:05.665
url: /archives/LeetCode-N-th Tribonacci Number
categories: 数据结构与算法
tags: 
---


Tribonacci序列Tn定义:

T0=0, T1=1, T2=1, n>=0时，Tn+3 = Tn + Tn+1 + Tn+2

限制条件是: 0<=n<=37, 32位整型。

我直接用C++撸了下面的代码，

```C++
#include <iostream>

using namespace std;

class Solution {
public:
    int tribonacci(int n) {
        if (n ==  2){
            return 1;
        } 
        if (n == 1){
            return 1;
        } 
        if (n == 0){
            return 0;
        } 
        int result = tribonacci(n-3) + tribonacci(n-2) + tribonacci(n-1);

        return result;
        
    }
};

// 下面是自己电脑编译测试用的代码
int main(){

	Solution sol;
	for (int i = 0; i <= 37; i++){
	    cout << "tribonacci " << i << " is " << sol.tribonacci(i) << endl;
	}

}
```

提交答案后，提示"The Limit Exceeded". 

在自己服务器测试时，也发现当n=31之后的速度下降的非常快。

```bash
# 编译
g++ -o trib trib.cpp
./trib
```

原因是当n越大时，直接使用递归会产生大量的重复计算，导致计算效率下降。为了解决该问题，需要维护一个已计算结果的字典，避免重复运算。

```c++
#include <iostream>

using namespace std;

class Solution {
public:
    int tribo[100];
    int tribonacci(int n) {
        if (n ==  2){
            return 1;
        } 
        if (n == 1){
            return 1;
        } 
        if (n == 0){
            return 0;
        } 
        
        if (tribo[n] > 0){
            return tribo[n]-1;
        }
        
        int result = tribonacci(n-3) + tribonacci(n-2) + tribonacci(n-1);
        tribo[n] = result + 1;

        return result;
        
    }
};

// 下面是自己电脑编译测试用的代码
int main(){

	Solution sol;
	for (int i = 0; i <= 37; i++){
	    cout << "tribonacci " << i << " is " << sol.tribonacci(i) << endl;
	}

}
```

新的代码中，我预先定义了一个大小为100的数组，主要是因为题目限制n的取值，如果n没有限制取值，我会考虑使用vector容器或者map词典。其次考虑到n=0时，结果为0，并且数组初始化值为0，为了能够通过数组值是否为0来判断是否已经计算了相应的tribonacci值，我对结果加上了一个pseducount。

最后这段代码运行速度是0ms.

如果要用哈希字典，代码如下

```C++
#include <map>
//一定要引入map头文件
class Solution {
public:
    //定义一个哈希字典
    map<int, int> num;

    int tribonacci(int n) {
        if(n < 2) return n;
        if(n == 2) return 1;
        //如果字典中有已经计算的结果,直接返回结果
        if(num.find(n) != num.end()) return num[n];
        int sum = tribonacci(n-1) + tribonacci(n-2) + tribonacci(n-3);
        num[n] = sum;
        return sum;
    }
};
```
