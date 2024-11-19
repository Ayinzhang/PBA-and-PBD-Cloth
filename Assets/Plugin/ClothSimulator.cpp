#include <cstdio>
#include <cmath>
#include <thread>
using namespace std;

struct constraint { int x, y; float z; constraint(int a = 0, int b = 0, float c = 0) : x(a), y(b), z(c) {}};

struct float3
{
    float x, y, z;
    float3(float a = 0, float b = 0, float c = 0) : x(a), y(b), z(c) {}
    float3 operator+(const float3& other) const { return float3(x + other.x, y + other.y, z + other.z); }
    float3 operator-(const float3& other) const { return float3(x - other.x, y - other.y, z - other.z); }
    friend float3 operator*(float other, const float3& vector) { return float3(vector.x * other, vector.y * other, vector.z * other); }
    float3& operator+=(const float3& other) { x += other.x; y += other.y; z += other.z; return *this; }
    float3& operator-=(const float3& other) { x -= other.x; y -= other.y; z -= other.z; return *this; }
    float Magnitude() const { return sqrt(x * x + y * y + z * z);}
    float3 Normalized() const { float mag = Magnitude(); return float3(x / mag, y / mag, z / mag); }
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

const int n = 35; bool init; float3 P1[n * n], P2[n * n], F[n * n]; constraint C[5 * n * n - 8 * n + 1]; long long CN[n * n];
thread th[n * n];

void Init(Mat4x4 mat, float3 g)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            P1[i * n + j] = mat * float3(2.0f * i / (n - 1) - 1, 2.0f * j / (n - 1) - 1, 0), F[i * n + j] = g;

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

void Transform(float3 Psrc[], float3 Pdst[], Mat4x4 mat, int l, int r)
{
    for (int i = l; i < r; i++) Pdst[i] = mat * Psrc[i];
}

void HandleConstaint(float3 P[], int k, int l, int r)
{
    for (int i = l; i < r; i++)
    {
        float3 p; long long c = CN[i];
        while (c > 0)
        {
            switch (c % 16)
            {
            case 1: p = P[i] - P[i + n]; F[i] -= k * (p.Magnitude() - 2.0f / (n - 1)) * p.Normalized(); break;
            case 2: p = P[i] - P[i - n]; F[i] -= k * (p.Magnitude() - 2.0f / (n - 1)) * p.Normalized(); break;
            case 3: p = P[i] - P[i + 1]; F[i] -= k * (p.Magnitude() - 2.0f / (n - 1)) * p.Normalized(); break;
            case 4: p = P[i] - P[i - 1]; F[i] -= k * (p.Magnitude() - 2.0f / (n - 1)) * p.Normalized(); break;
            case 5: p = P[i] - P[i + n + 1]; F[i] -= k * (p.Magnitude() - 2.828f / (n - 1)) * p.Normalized(); break;
            case 6: p = P[i] - P[i - n - 1]; F[i] -= k * (p.Magnitude() - 2.828f / (n - 1)) * p.Normalized(); break;
            case 7: p = P[i] - P[i + n - 1]; F[i] -= k * (p.Magnitude() - 2.828f / (n - 1)) * p.Normalized(); break;
            case 8: p = P[i] - P[i - n + 1]; F[i] -= k * (p.Magnitude() - 2.828f / (n - 1)) * p.Normalized(); break;
            case 9: p = P[i] - P[i + 2 * n]; F[i] -= k * (p.Magnitude() - 4.0f / (n - 1)) * p.Normalized(); break;
            case 10: p = P[i] - P[i - 2 * n]; F[i] -= k * (p.Magnitude() - 4.0f / (n - 1)) * p.Normalized(); break;
            case 11: p = P[i] - P[i + 2]; F[i] -= k * (p.Magnitude() - 4.0f / (n - 1)) * p.Normalized(); break;
            case 12: p = P[i] - P[i - 2]; F[i] -= k * (p.Magnitude() - 4.0f / (n - 1)) * p.Normalized(); break;
            }
            c >>= 4;
        }
    }
}

void UpdatePosition(float3 P1[], float3 P2[], float3 g, float d, float t2, int l, int r)
{
    for (int i = l; i < r; i++)
    {
        if (i != n - 1 && i != n * n - 1) P2[i] += d * (P2[i] - P1[i]) + t2 * F[i];
        P1[i] = P2[i]; F[i] = g;
    }
}

extern "C"  _declspec(dllexport) void Semi_Implict_Cpp(float3 P[], float3 P1[], float3 P2[],
    Mat4x4 mat, Mat4x4 imat, float3 g, int k, int iter, float dt, float damping)
{
    if(!init) Init(mat, g), init = true;

    for (int i = 0; i < n * n; i++) P2[i] = mat * P[i];

    float t = dt / iter, t2 = t * t, d = pow(damping, 1.0f / iter);
    for (int i = 0; i < iter; i++)
    {
        // Constrant
        for (int j = 0; j < 5 * n * n - 8 * n + 1; j++)
        {
            int x = C[j].x, y = C[j].y; float z = C[j].z;
            float3 f = k * ((P2[x] - P2[y]).Magnitude() - z) * (P2[x] - P2[y]).Normalized();
            F[x] -= f; F[y] += f;
        }

        // Gravity & Update & Calculate Local Position
        for (int j = 0; j < n * n; j++)
        {
            if (j != n - 1 && j != n * n - 1) P2[j] += d * (P2[j] - P1[j]) + t2 * F[j];
            P1[j] = P2[j]; F[j] = g;
        }
    }

    for (int i = 0; i < n * n; i++) P[i] = imat * P2[i];
}

extern "C"  _declspec(dllexport) void Semi_Implict_CppThread(float3 P[], float3 P1[], float3 P2[],
    Mat4x4 mat, Mat4x4 imat, float3 g, int k, int iter, int thr, float dt, float damping)
{
    if (!init) Init(mat, g), init = true;

    for (int i = 0; i < thr; i++) th[i] = thread(Transform, P, ref(P2), mat, n * n / thr * i, i == thr - 1 ? n * n + 1: n * n / thr * (i + 1));
    for (int i = 0; i < thr; i++) th[i].join();

    float t = dt / iter, t2 = t * t, d = pow(damping, 1.0f / iter);
    for (int i = 0; i < iter; i++)
    {
        // Constrant
        for (int i = 0; i < thr; i++) th[i] = thread(HandleConstaint, P2, k, n * n / thr * i, i == 15 ? n * n + 1 : n * n / thr * (i + 1));
        for (int i = 0; i < thr; i++) th[i].join();

        // Gravity & Update & Calculate Local Position
        for (int i = 0; i < thr; i++) th[i] = thread(UpdatePosition, ref(P1), ref(P2), g, d, t2, n * n / thr * i, i == 15 ? n * n + 1 : n * n / thr * (i + 1));
        for (int i = 0; i < thr; i++) th[i].join();
    }

    for (int i = 0; i < thr; i++) th[i] = thread(Transform, P2, ref(P), imat, n * n / thr * i, i == 15 ? n * n + 1 : n * n / thr * (i + 1));
    for (int i = 0; i < thr; i++) th[i].join();

}