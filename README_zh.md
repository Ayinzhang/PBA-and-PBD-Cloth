# PBA和PBD布料

[English Version](README.md)

*学校座谈会需要选取论文演讲，加之之前想整PBD，XPBD相关但一直卡在读论文阶段，就开了这个项目。旨在对比运行环境(Unity，Unity Job+Burst，C++, C++Thread)以及模拟方法(半隐式欧拉，隐式欧拉，PBD，XPBD)对模拟性能的影响。（以下例子参数为21x21布料网格，2000劲度系数，32迭代，64线程，GTX1650笔记本。并行均为无锁多线程）*

|           | Unity | Unity Job + Burst | C++  | C++ Parallelize |
| --------- | ----- | ----------------- | ---- | --------------- |
| Time (ms) | 5.6   | 0.47              | 0.56 | 0.76            |

*C++ Thread的表现有些出乎意料，感觉可能是Thread过于耗了？不过C++比Unity的性能倒好了不少，感觉之后遇到想优化不适合Job + Burst并行化的代码都可以试试C++*

|                  | Semi-Implicit Euler | Implicit Euler | PBD  | XPBD |
| ---------------- | ------------------- | -------------- | ---- | ---- |
| Main Thread (ms) | 5.6                 | 6.0            | 11.2 | 8.7  |
| Job + Burst (ms) | 0.47                | 0.58           | 0.48 | 0.68 |

*隐式的炸了不计入考量，主线程下XPBD优于PBD应该是我在写PBD的时候用的雅可比而非高斯塞德尔，但看并行化的就行了。比较不甘的是Unity自己的Cloth组件模拟这些只需要0.2ms，而我的PBD在保持稳定的情况下只能到0.4ms，小作坊下猛料也没下过官方，之后有空再研究下如何优化。*

PBD Cloth:

![img](https://pica.zhimg.com/80/v2-f1a0758cdaa72daa7836aec27a6ac8df_720w.gif?source=d16d100b)

XPBD Cloth:

![img](https://picx.zhimg.com/80/v2-fb2d697cbb11746a20b37886ba0fc903_720w.gif?source=d16d100b)
