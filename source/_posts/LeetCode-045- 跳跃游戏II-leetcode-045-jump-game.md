---
title: LeetCode-045- 跳跃游戏II
date: 2020-03-10 14:46:50.473
updated: 2020-03-10 14:47:07.849
url: /archives/leetcode-045-jump-game
categories: 数据结构与算法
tags: 数据结构 | C/C++
---

# 跳跃游戏II

> [https://leetcode-cn.com/problems/jump-game-ii/](https://leetcode-cn.com/problems/jump-game-ii/)

相对与之前的跳跃游戏，这道题目保证能够抵达终点，但是要求你输出最短路径。

DP三部曲，定义子问题，定义状态数组，和状态转移方程。

- 子问题，到i个位置的最短步数，等于0..i-1中能够到i的最小步数+1
- 状态数组: dp[i]
- 状态转移方程: dp[i] = min(dp[i],dp[j]+1), j =0..i

下面就是最初版的代码

```cpp
// 90 / 92 个通过测试用例
class Solution1 {
public:
    int jump(vector<int>& nums) {
        int dp[nums.size()];
        for (int i = 0; i < nums.size(); i++) dp[i] = INT_MAX;
        dp[0] = 0; //初始化第一个位置
        for (int i = 1 ;i < nums.size() ; i++){
            for(int j = 0; j < i ; j++){
                if (nums[j] + j >= i){
                    dp[i] = min(dp[i],dp[j]+1);
                }
            }
        }
        return dp[nums.size()-1];
    }
};
```

虽然代码正确，也能得到正确的结果，但是超时了，卡在了一个全为1的输入中。

**优化下代码**：外层依旧遍历整个数组，内层循环则是遍历当前数组能够抵达的位置， 更改其最小步数。

也就是状态转移方程变成了dp[i+j] = min(dp[i+j], dp[i]+1);

这就保证了对于全为1的数组，内层循环只会走一步，就避免了之前的情况。

```cpp
class Solution2 {
public:
    int jump(vector<int>& nums) {
        int dp[nums.size()];
        int num_len = nums.size() ;
        for (int i = 0; i < nums.size(); i++) dp[i] = INT_MAX;
        dp[0] = 0; //初始化第一个位置

        for (int i = 0 ;i < nums.size() ; i++){
            for(int j = 1; j <= nums[i] ; j++){
                int next = min(i+j, num_len-1); //不能越界
                dp[next] = min(dp[next],dp[i]+1);
            }
        }
        return dp[nums.size()-1];
    }
};
```

当我满怀信心运行代码的时候，结果还是超时了。

原来输入中，还有一个单调递减的情况，导致最终内层循环次数分别为n,n-1,n-2...1位置，导致许多不必要的运算。

**继续优化**：我们新增一个变量，farthest，表示每次能够抵达的最远距离。如果不能保证当前位置能够比farthest走的更远的时候, 就说明当前位置不够好。

同时，在内层循环之前，我们先判断当前位置是否能够直接抵达终点，如果可以就直接返回结果。

```cpp
class Solution {
public:
    int jump(vector<int>& nums) {
        int num_len = nums.size() ;
        if ( num_len <= 1) return 0;
        int farthest = 0;
        int dp[num_len];
        for (int i = 0; i < num_len ; i++) dp[i] = INT_MAX - 1;
        dp[0] = 0; //初始化第一个位置

        for (int i = 0 ;i < num_len ; i++){
            if ((i+nums[i]+1)>= num_len) return dp[i]+1;
            for(int j = 1; j <= nums[i] && nums[i] + i > farthest; j++){
                int next = i + j; //抵达的最后位置
                dp[next] = min(dp[next],dp[i]+1);
            }
            farthest = i + nums[i];
        }
        return dp[num_len-1];
    }
};
```

**继续优化**:  我们内层循环每次也不是需要从**当前的下一个位置**开始更新，而是farthest-i开始，因为next=i+j = i+farthest-i+j = farthest, 即我们从**之前所能抵达最远位置**开始更新，一直更新到**当前所能抵达最远位置**。

```cpp
class Solution {
public:
    int jump(vector<int>& nums) {
        int num_len = nums.size() ;
        if ( num_len <= 1) return 0;
        int farthest = 0;
        int dp[num_len];
        for (int i = 0; i < num_len ; i++) dp[i] = INT_MAX - 1;
        dp[0] = 0; //初始化第一个位置

        for (int i = 0 ;i < num_len ; i++){
            if ((i+nums[i]+1)>= num_len) return dp[i]+1;
            for(int j = farthest - i; j <= nums[i] && nums[i] + i > farthest; j++){
                int next = i + j; //抵达的最后位置
                dp[next] = min(dp[next],dp[i]+1);
            }
            farthest = i + nums[i];
        }
        return dp[num_len-1];
    }
};
```

最后的代码提交之后运行速度是4ms。

同样是两层循环，通过剪枝能够大大提高运算速度。