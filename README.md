# PBA and PBD Cloth

[中文版](README_zh.md)

*The school symposium requires the paper presentations. In addition, I wanted to study PBD and XPBD before but was stuck at the paper reading stage, so I started this project. It aims to compare the impact of the operating environment (Unity, Unity Job+Burst, C++, C++Thread) and simulation methods (Semi-implicit Euler, Implicit Euler, PBD, XPBD) on performance. (The following example parameters are 35x35 cloth grid, 2000 stiffness coefficient, 32 iterations, 64 threads, GTX1650 laptop. Parallelism is all lock-free multi-threading)*

|           | Unity | Unity Job + Burst | C++   | C++ Thread |
| --------- | ----- | ----------------- | ----- | ---------- |
| Time (ms) | 36.9  | 1.6               | 12.2  | 38.8       |

*The performance of C++ Thread is a bit unexpected. I think it might be because Thread too costly? However, the performance of C++ is much better than Unity. I think I can try C++ when I want to optimize the code that is not suitable for Job + Burst parallelization the next time.*

|                  | Semi-Implicit Euler | Implicit Euler | PBD  | XPBD |
| ---------------- | ------------------- | -------------- | ---- | ---- |
| Main Thread (ms) | 36.9                | N/A            | 70.9 | 52.2 |
| Job + Burst (ms) | 1.6                 | N/A            | 2.1  | 2.2  |

*Implicit crashes and not be taken into account. XPBD is better than PBD in the main thread may because I used Jacobi instead of Gauss-Seidel when writing PBD, but just look at the parallelization. What disappoint me is that Unity's Cloth component only needs 0.2ms to simulate these, while my XPBD can only reach 0.4ms while maintaining stability. I will study how to optimize it later when I have time.* 
