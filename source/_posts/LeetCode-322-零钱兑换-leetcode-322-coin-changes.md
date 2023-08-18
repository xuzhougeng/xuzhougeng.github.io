---
title: LeetCode-322-零钱兑换
date: 2020-03-09 18:13:26.216
updated: 2020-03-09 18:29:42.993
url: /archives/leetcode-322-coin-changes
categories: 数据结构与算法
tags: C/C++
---

# 硬币兑换

> 来源LeetCode, 题目地址<<https://leetcode-cn.com/problems/coin-change/>>

给定不同面额的硬币 coins 和一个总金额 amount。编写一个函数来计算可以凑成总金额所需的最少的硬币个数。如果没有任何一种硬币组合能组成总金额，返回 -1。

## 暴力递归（回溯）

每一种硬币都有选择和不选两种情况，因此最简单的方法就是穷举所有可能性，从中找到最小的选择，递归树如下:

![状态树](https://halo-1252249331.cos.ap-shanghai.myqcloud.com/upload/2020/3/image-e12cb3fbb81b4bbea59758f259b7eb2b.png)

每次节点的子节点都会出现3个分支，因此时间复杂度是指数级。

尽管明知这个代码是不能通过LeetCode，我还是把代码给写完了。

```cpp
class Solution1 {
public:
    int MIN_COMB = INT_MAX;
    int coinChange(vector<int>& coins, int amount) {
        DFS(0, amount, coins);
        return MIN_COMB;
    }

    void DFS(int depth, int amount, vector<int>& coins){
        if (amount == 0){
            MIN_COMB = min(depth, MIN_COMB);
            return ;
        }
        for (auto coin : coins){
            if (amount - coin < 0 ) continue; //减枝
            DFS(depth+1, amount - coin, coins);
            
        }
        return ;
    }
};
```

对于一些比较正常的数据，这个代码都能正常运行，但是对于LeetCode中的`[186,419,83,408] 6249`，标准答案是20，也就是至少有20层，因此计算量至少是`3^20`，但是递归深度可以达到75层，计算量就是一个天文数字了。

## 递归的广度优先搜索

递归树其实是一种深度优先的搜索方法，我们可以采用广度优先搜索方法，对递归树按层进行遍历，第一次遇到amount=0的时候，也就是最小的遍历层数。

```cpp
class Solution2 {
public:
    int coinChange(vector<int>& coins, int amount) {

        queue<pair<int,int>> q; //level and amount
        q.push({0, amount});

        while ( !q.empty()) {

            int level = q.front().first;
            int surplus = q.front().second;
            q.pop();

            for (auto c : coins){
                if ( surplus - c < 0 ) continue;
                if ( surplus - c == 0 ) return level + 1;
                q.push({level+1, surplus-c});
            }
        }
        return -1;


    }
};
```

尽管看起来， BFS或许比DFS能够更快的找到最终的解（不需要遍历所有的状态），但是由于递归树本身就是指数增长，因此只要层数过大，时间就会惊人的增长。

你可以保持`[1,2,5]`不变，然后依次测试，10,20,50,80,100. 你会发现前4个解决速度都比较快，但是到了100，时间就上天了。

## 贪心算法

暴力递归的问题在于，节点的增长是指数级别的。如果我们每次都选择当前的最优解，也就是对于11，从1，2，5中找交换的硬币的时候，都以5，2，1这种顺序进行选择，就避免了指数级别的增长。

```cpp
class Solution {
public:
    int coinChange(vector<int>& coins, int amount) {
        sort(coins.begin(), coins.end(), std::greater<>());
        return DFS(0, amount, coins);
    }

    int DFS(int level, int amount, vector<int>& coins){
        if (amount <= 0){
            return amount == 0 ? level : -1;
        }
        for (auto coin : coins){
            int res = DFS(level+1, amount - coin, coins);
            if (res > 0 ) return res;
        }
        return -1;
    }
};
```

但是贪心算法有一个弊端，你得要证明能够从局部最优一直推导出全局最优解，才能保证你的结果是正确的。如果无法保证，那么贪心算法不一定保证最终结果就是最优解，所以这里贪心算法也是无法通过。

## 递归+记忆化

在我们之前的递归树中，无论采用何种遍历方式（DFS或BFS），都无法逃脱时间指数级别增长的命运。

不过当我们仔细观察递归树的时候，我们会发现一些节点是被重复计算的。如果能避免这些重复计算，那么时间复杂度就会降低到O(n^2)

在原来递归的代码上加上记忆化，我写了很久，主要是不知道在递归的时候，如何记录当前情况下，如何挑选最优子问题。

先放上正确的答案

```cpp
class Solution {
public:
    unordered_map<int, int> dict;
    int coinChange(vector<int>& coins, int amount) {
        if (amount < 1) return 0;
        return DFS(coins, amount);
    }
    int DFS(vector<int>& coins, int amount){
        if (amount < 0 ) return -1;
        if (amount == 0 ) return 0;
        if (dict.find(amount) != dict.end()) return dict[amount];

        int min = INT_MAX; //足够大的值
        for (int coin : coins){
            int res = DFS(coins, amount-coin);
            if (res >= 0 && res < min){
                min = res + 1;
            }
        }
        dict[amount] = (min == INT_MAX ? -1 : min);
        return dict[amount];
    }
};
```

递归返回的从最后一层到当前层的所需的步数。而之前的递归程代码则是每次记录当前的深度，最终拿深度和全局最小值进行比较。

下面则是我的错误示范。我的代码主题也是想算出最后位置到当前位置的经历的步数，但是我想用最后一层的深度减去当前的深度。结果每一个amount对应都是1.

```cpp
// 递归+记忆化
class Solution4 {
public:
    unordered_map<int, int> dict;
    int coinChange(vector<int>& coins, int amount) {
    
        return DFS(0, amount, coins);
    }

    int DFS(int depth, int amount, vector<int>& coins){
        if (amount == 0){
            return depth;
        }

        if (dict.find(amount) != dict.end() ) {
            return dict[amount];
        }

        int local_min = INT_MAX;

        for (auto coin : coins){
            if (amount - coin < 0 ) continue; 
            local_min = min(local_min, DFS(depth+1, amount - coin, coins) -depth);
            
        }
        dict[amount] = (local_min == INT_MAX ? -1 :   local_min);
        return dict[amount];
    }
};
```

## 自底向上的动态规划

事实上递归加记忆化就是一种动态规划，只不过递归是一种自顶向下的策略。并且有了之前递归加记忆化的经验，我们写递推就变得简单了。

第一步: 定义子问题。我们求解组成金额为n的最少硬币数，也就是求解一系列 n-k (k=coins) 子问题中的最小选择加1。以amount=11，coins=1,2,5为例。求解amount=11的最少硬币，也就是在10,9,6这三种选择中挑选其中所需硬币最小的子问题，然后加1.

第二步：定义DP数组. f(n) = min{f(n-k), k = 1,2,5} + 1

第三步:  定义DP方程: dp[i] = min(dp[i], dp[i-coins[j] ] )+ 1

最终代码如下

```cpp
class Solution {
public:
    int coinChange(vector<int>& coins, int amount) {
        int MAX = amount + 1;// 只要保证比amount大即可, 因为后续要和当前最小的选择比较。
        vector<int> dp(amount+1, MAX); //DP数组
        dp[0] = 0;
        for (int i = 1; i <= amount; i++){
            for (int j = 0; j < coins.size(); j++){
                if (coins[j] <= i){ //举例, i = 2, 只能考虑coin=1,2, 排除5
                    dp[i] = min(dp[i], dp[i-coins[j]] + 1); 
                }
            }
        }
        return dp[amount] > amount ? -1 : dp[amount];
    }
};
```

代码中的细节，

- 初始化的数组大小为amount+1,  不能是amount, 否则只能表示0到amount-1
- 初始化的数组存放的值，要比最坏情况下的最少硬币数大，例如amount=10, 硬币只有11，那么数组内容就会是原本情况，最终的结果是比amount大，说明无解。因此，MAX = INT\_MAX-1 也是可以的。

最终来看，自底向上的动态递归的代码反而是最简洁的。