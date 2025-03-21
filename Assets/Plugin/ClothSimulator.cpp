#include <omp.h>
#include <cstdio>
#include <cmath>
using namespace std;

struct constraint { int x, y; float z; constraint(int a = 0, int b = 0, float c = 0) : x(a), y(b), z(c) {} };

struct float3
{
    float x, y, z;
    float3(float a = 0, float b = 0, float c = 0) : x(a), y(b), z(c) {}
    float3 operator+(const float3& other) const { return float3(x + other.x, y + other.y, z + other.z); }
    float3 operator-(const float3& other) const { return float3(x - other.x, y - other.y, z - other.z); }
    friend float3 operator*(float other, const float3& vector) { return float3(vector.x * other, vector.y * other, vector.z * other); }
    friend float3 operator/(const float3& vector, float other) { return float3(vector.x / other, vector.y / other, vector.z / other); }
    float3& operator+=(const float3& other) { x += other.x; y += other.y; z += other.z; return *this; }
    float3& operator-=(const float3& other) { x -= other.x; y -= other.y; z -= other.z; return *this; }
    float Magnitude() const { return sqrt(x * x + y * y + z * z); }
};

struct Mat4x4
{
    float thr[16];

    float3 operator*(const float3& vec) const // Matrices in Unity are in column-major order
    {
        float x = thr[0] * vec.x + thr[4] * vec.y + thr[8] * vec.z + thr[12];
        float y = thr[1] * vec.x + thr[5] * vec.y + thr[9] * vec.z + thr[13];
        float z = thr[2] * vec.x + thr[6] * vec.y + thr[10] * vec.z + thr[14];
        float w = thr[3] * vec.x + thr[7] * vec.y + thr[11] * vec.z + thr[15];

        if (w != 0) { x /= w; y /= w; z /= w; }
        return float3(x, y, z);
    }
};

const int n = 21; bool init; float3 F[n * n]; constraint C[5 * n * n - 8 * n + 1]; long long CN[n * n]; thread th[n * n];

void Init()
{
    int cnt = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            // Stretch
            if (i < n - 1)
            {
                C[cnt++] = constraint(i * n + j, (i + 1) * n + j, 2.0f / (n - 1));
                CN[i * n + j] = (CN[i * n + j] << 4) + 1; CN[(i + 1) * n + j] = (CN[(i + 1) * n + j] << 4) + 2;
            }
            if (j < n - 1)
            {
                C[cnt++] = constraint(i * n + j, i * n + j + 1, 2.0f / (n - 1));
                CN[i * n + j] = (CN[i * n + j] << 4) + 3; CN[i * n + j + 1] = (CN[i * n + j + 1] << 4) + 4;
            }

            // Shear
            if ((i + j) % 2 == 0)
            {
                if (i < n - 1 && j < n - 1)
                {
                    C[cnt++] = constraint(i * n + j, (i + 1) * n + j + 1, 2.828f / (n - 1));
                    CN[i * n + j] = (CN[i * n + j] << 4) + 5; CN[(i + 1) * n + j + 1] = (CN[(i + 1) * n + j + 1] << 4) + 6;
                }
                if (i < n - 1 && j > 0)
                {
                    C[cnt++] = constraint(i * n + j, (i + 1) * n + j - 1, 2.828f / (n - 1));
                    CN[i * n + j] = (CN[i * n + j] << 4) + 7; CN[(i + 1) * n + j - 1] = (CN[(i + 1) * n + j - 1] << 4) + 8;
                }
            }

            // Bend
            if (i < n - 2)
            {
                C[cnt++] = constraint(i * n + j, (i + 2) * n + j, 4.0f / (n - 1));
                CN[i * n + j] = (CN[i * n + j] << 4) + 9; CN[(i + 2) * n + j] = (CN[(i + 2) * n + j] << 4) + 10;
            }
            if (j < n - 2)
            {
                C[cnt++] = constraint(i * n + j, i * n + j + 2, 4.0f / (n - 1));
                CN[i * n + j] = (CN[i * n + j] << 4) + 11; CN[i * n + j + 2] = (CN[i * n + j + 2] << 4) + 12;
            }
        }
}

void TransformIn(float3 Psrc[], float3 Pdst[], float3 Plast[], float3 V[], Mat4x4 mat, float dt, int l, int r)
{
    for (int i = l; i < r; i++) Pdst[i] = mat * Psrc[i], V[i] = (Pdst[i] - Plast[i]) / dt;
}

void TransformOut(float3 Psrc[], float3 Pdst[], Mat4x4 mat, int l, int r)
{
    for (int i = l; i < r; i++) Pdst[i] = mat * Psrc[i];
}

void HandleConstaint(float3 P[], long long CN[], int k, int l, int r)
{
    for (int i = l; i < r; i++)
    {
        float3 p; long long c = CN[i];
        while (c > 0)
        {
            switch (c % 16)
            {
            case 1: p = P[i] - P[i + n]; F[i] -= k * (1 - 2.0f / (n - 1) / p.Magnitude()) * p; break;
            case 2: p = P[i] - P[i - n]; F[i] -= k * (1 - 2.0f / (n - 1) / p.Magnitude()) * p; break;
            case 3: p = P[i] - P[i + 1]; F[i] -= k * (1 - 2.0f / (n - 1) / p.Magnitude()) * p; break;
            case 4: p = P[i] - P[i - 1]; F[i] -= k * (1 - 2.0f / (n - 1) / p.Magnitude()) * p; break;
            case 5: p = P[i] - P[i + n + 1]; F[i] -= k * (1 - 2.828f / (n - 1) / p.Magnitude()) * p; break;
            case 6: p = P[i] - P[i - n - 1]; F[i] -= k * (1 - 2.828f / (n - 1) / p.Magnitude()) * p; break;
            case 7: p = P[i] - P[i + n - 1]; F[i] -= k * (1 - 2.828f / (n - 1) / p.Magnitude()) * p; break;
            case 8: p = P[i] - P[i - n + 1]; F[i] -= k * (1 - 2.828f / (n - 1) / p.Magnitude()) * p; break;
            case 9: p = P[i] - P[i + 2 * n]; F[i] -= k * (1 - 4.0f / (n - 1) / p.Magnitude()) * p; break;
            case 10: p = P[i] - P[i - 2 * n]; F[i] -= k * (1 - 4.0f / (n - 1) / p.Magnitude()) * p; break;
            case 11: p = P[i] - P[i + 2]; F[i] -= k * (1 - 4.0f / (n - 1) / p.Magnitude()) * p; break;
            case 12: p = P[i] - P[i - 2]; F[i] -= k * (1 - 4.0f / (n - 1) / p.Magnitude()) * p; break;
            }
            c >>= 4;
        }
    }
}

void UpdatePosition(float3 P1[], float3 P2[], float3 V[], float3 g, float d, float idt, int l, int r)
{
    for (int i = l; i < r; i++)
    {
        V[i] = d * V[i] + idt * F[i];
        if (i != n - 1 && i != n * n - 1) P2[i] += idt * V[i];
        P1[i] = P2[i]; F[i] = g;
    }
}

extern "C"  _declspec(dllexport) void Semi_Implict_Cpp(float3 P[], float3 P1[], float3 P2[], float3 V[],
    Mat4x4 mat, Mat4x4 imat, float3 g, int k, int iter, float dt, float damping)
{
    if (!init) Init(), init = true;

    for (int i = 0; i < n * n; i++) P2[i] = mat * P[i], V[i] += (P2[i] - P1[i]) / dt;

    float idt = dt / iter, d = pow(damping, 1.0f / iter);
    for (int i = 0; i < iter; i++)
    {
        // Constrant
        for (int j = 0; j < 5 * n * n - 8 * n + 1; j++)
        {
            int x = C[j].x, y = C[j].y; float z = C[j].z;
            float3 f = k * (1 - z / (P2[x] - P2[y]).Magnitude()) * (P2[x] - P2[y]);
            F[x] -= f; F[y] += f;
        }

        // Gravity & Update & Calculate Local Position
        for (int j = 0; j < n * n; j++)
        {
            V[j] = d * V[j] + idt * F[j];
            if (j != n - 1 && j != n * n - 1) P2[j] += idt * V[j];
            P1[j] = P2[j]; F[j] = g;
        }
    }

    for (int i = 0; i < n * n; i++) P[i] = imat * P2[i];
}

extern "C" _declspec(dllexport) void Semi_Implict_CppThread(float3 P[], float3 P1[], float3 P2[], float3 V[],
    Mat4x4 mat, Mat4x4 imat, float3 g, int k, int iter, int thr, float dt, float damping)
{
    if (!init) Init(), init = true;

#pragma omp parallel num_threads(thr)
    {
        int tid = omp_get_thread_num();
        int nthr = omp_get_num_threads();
        int total = n * n;
        int chunk = total / nthr;
        int start = tid * chunk;
        int end = (tid == nthr - 1) ? total : start + chunk;
        TransformIn(P, P2, P1, V, mat, dt, start, end);
    }

    float idt = dt / iter, d = pow(damping, 1.0f / iter);
    for (int i = 0; i < iter; i++)
    {
#pragma omp parallel num_threads(thr)
        {
            int tid = omp_get_thread_num();
            int nthr = omp_get_num_threads();
            int total = n * n;
            int chunk = total / nthr;
            int start = tid * chunk;
            int end = (tid == nthr - 1) ? total : start + chunk;
            HandleConstaint(P2, CN, k, start, end);
        }

#pragma omp parallel num_threads(thr)
        {
            int tid = omp_get_thread_num();
            int nthr = omp_get_num_threads();
            int total = n * n;
            int chunk = total / nthr;
            int start = tid * chunk;
            int end = (tid == nthr - 1) ? total : start + chunk;
            UpdatePosition(P1, P2, V, g, d, idt, start, end);
        }
    }

#pragma omp parallel num_threads(thr)
    {
        int tid = omp_get_thread_num();
        int nthr = omp_get_num_threads();
        int total = n * n;
        int chunk = total / nthr;
        int start = tid * chunk;
        int end = (tid == nthr - 1) ? total : start + chunk;
        TransformOut(P2, P, imat, start, end);
    }
}