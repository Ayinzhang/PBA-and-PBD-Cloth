# PBA and PBD Cloth

[中文版](README_zh.md)

*The school symposium requires the paper presentations. In addition, I wanted to study PBD and XPBD before but was stuck at the paper reading stage, so I started this project. It aims to compare the impact of the operating environment (Unity, Unity Job+Burst, C++, C++Thread) and simulation methods (Semi-implicit Euler, Implicit Euler, PBD, XPBD) on performance. (The following example parameters are 21x21 cloth grid, 2000 stiffness coefficient, 32 iterations, 64 threads, GTX1650 laptop. Parallelism is all lock-free multi-threading)*

|           | Unity | Unity Job + Burst | C++  | C++ Parallelize |
| --------- | ----- | ----------------- | ---- | --------------- |
| Time (ms) | 5.6   | 0.47              | 0.56 | 0.76            |

*The performance of C++ Thread is a bit unexpected. I think it might be because Thread too costly? However, the performance of C++ is much better than Unity. I think I can try C++ when I want to optimize the code that is not suitable for Job + Burst parallelization the next time.*

|                  | Semi-Implicit Euler | Implicit Euler | PBD  | XPBD |
| ---------------- | ------------------- | -------------- | ---- | ---- |
| Main Thread (ms) | 5.6                 | 6.0            | 11.2 | 8.7  |
| Job + Burst (ms) | 0.47                | 0.58           | 0.48 | 0.68 |

*Implicit crashes and not be taken into account. XPBD is better than PBD in the main thread may because I used Jacobi instead of Gauss-Seidel when writing PBD, but just look at the parallelization. What disappoint me is that Unity's Cloth component only needs 0.2ms to simulate these, while my PBD can only reach 0.4ms while maintaining stability. I will study how to optimize it later when I have time.* 

PBD Cloth:

![img](https://pica.zhimg.com/80/v2-f1a0758cdaa72daa7836aec27a6ac8df_720w.gif?source=d16d100b)

XPBD Cloth:

![img](https://picx.zhimg.com/80/v2-fb2d697cbb11746a20b37886ba0fc903_720w.gif?source=d16d100b)
