using System;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.Jobs;
using Unity.Jobs; // Define¡¾IJob¡¿¡¾IJobParallelFor¡¿
using Unity.Burst; // Define¡¾BurstCompile¡¿
using Unity.Collections; // Define¡¾NativeArray¡¿
using Unity.Mathematics; // SIMD Math Library
using Unity.Collections.LowLevel.Unsafe; // Quick Memory Copy
using System.Runtime.InteropServices;

public class ClothSimulator : MonoBehaviour
{
    public enum Type { Semi_Implict, Semi_Implict_Burst, Semi_Implict_Cpp, Semi_Implict_CppThread, Implict, PBD, XPBD}
    public Type type = Type.Semi_Implict; public int k = 5000, iter = 7, thread = 2; public float dt = 1f / 30, damping = 0.95f;

    Mesh mesh; Text[] texts; Slider[] sliders; Vector3 g = new Vector3(0, -9.8f, 0); float t; const int n = 35;
    Vector3[] P = new Vector3[n * n], P1 = new Vector3[n * n], P2 = new Vector3[n * n], F = new Vector3[n * n];
    Vector2[] UV = new Vector2[n * n]; int[] T = new int[(n - 1) * (n - 1) * 6];
    struct constraint { public int x, y; public float z; } constraint[] C = new constraint[5 * n * n - 8 * n + 1];

    JobHandle handle; NativeArray<float3> FN = new NativeArray<float3>(n * n, Allocator.Persistent), 
    P1N = new NativeArray<float3>(n * n, Allocator.Persistent), P2N = new NativeArray<float3>(n * n, Allocator.Persistent),
    P3N = new NativeArray<float3>(n * n, Allocator.Persistent); NativeArray<long> CN = new NativeArray<long>(5 * n * n - 8 * n + 1, Allocator.Persistent);

    [DllImport("ClothSimulator", CallingConvention = CallingConvention.Cdecl)]
    static extern void Semi_Implict_Cpp( [In, Out] Vector3[] P, [In, Out] Vector3[] P1, [In, Out] Vector3[] P2,
        Matrix4x4 mat, Matrix4x4 imat, Vector3 g, int k, int iter, float dt, float damping);

    [DllImport("ClothSimulator", CallingConvention = CallingConvention.Cdecl)]
    static extern void Semi_Implict_CppThread([In, Out] Vector3[] P, [In, Out] Vector3[] P1, [In, Out] Vector3[] P2,
        Matrix4x4 mat, Matrix4x4 imat, Vector3 g, int k, int iter, int thr, float dt, float damping);

    void Start()
    {
        mesh = GetComponent<MeshFilter>().mesh;
        texts = GameObject.Find("Canvas").GetComponentsInChildren<Text>();
        sliders = GameObject.Find("Canvas").GetComponentsInChildren<Slider>();

        // Resize the mesh (Blender type)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                P[i * n + j] = new Vector3(2f * i / (n - 1) - 1, 2f * j / (n - 1) - 1, 0);
                P1[i * n + j] = P1N[i * n + j] = transform.TransformPoint(P[i * n + j]);
                F[i * n + j] = FN[i * n + j] = g;
                UV[i * n + j] = new Vector2(1 - i / (n - 1f), j / (n - 1f));
            }
        int cnt = 0;
        for (int i = 0; i < n - 1; i++)
            for (int j = 0; j < n - 1; j++)
            {
                T[cnt * 6 + 0] = i * n + j;
                T[cnt * 6 + 1] = (i + 1) * n + j + 1;
                T[cnt * 6 + 2] = i * n + j + 1;
                T[cnt * 6 + 3] = i * n + j;
                T[cnt * 6 + 4] = (i + 1) * n + j;
                T[cnt++ * 6 + 5] = (i + 1) * n + j + 1;
            }
        mesh.Clear(); mesh.vertices = P; mesh.triangles = T;
        mesh.uv = UV; mesh.RecalculateNormals();

        // Create Constraint
        if (C[0].z > 0) return;
        cnt = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                // Stretch
                if (i < n - 1)
                {
                    C[cnt++] = new constraint { x = i * n + j, y = (i + 1) * n + j, z = 2f / (n - 1) };
                    CN[i * n + j] = (CN[i * n + j] << 4) + 1; CN[(i + 1) * n + j] = (CN[(i + 1) * n + j] << 4) + 2;
                }
                if (j < n - 1)
                {
                    C[cnt++] = new constraint { x = i * n + j, y = i * n + j + 1, z = 2f / (n - 1) };
                    CN[i * n + j] = (CN[i * n + j] << 4) + 3; CN[i * n + j + 1] = (CN[i * n + j + 1] << 4) + 4;
                }

                // Shear
                if ((i + j) % 2 == 0)
                {
                    if (i < n - 1 && j < n - 1)
                    {
                        C[cnt++] = new constraint { x = i * n + j, y = (i + 1) * n + j + 1, z = 2.828f / (n - 1) };
                        CN[i * n + j] = (CN[i * n + j] << 4) + 5; CN[(i + 1) * n + j + 1] = (CN[(i + 1) * n + j + 1] << 4) + 6;
                    }
                    if (i < n - 1 && j > 0)
                    {
                        C[cnt++] = new constraint { x = i * n + j, y = (i + 1) * n + j - 1, z = 2.828f / (n - 1) };
                        CN[i * n + j] = (CN[i * n + j] << 4) + 7; CN[(i + 1) * n + j - 1] = (CN[(i + 1) * n + j - 1] << 4) + 8;
                    }
                }

                // Bend
                if (i < n - 2)
                {
                    C[cnt++] = new constraint { x = i * n + j, y = (i + 2) * n + j, z = 4f / (n - 1) };
                    CN[i * n + j] = (CN[i * n + j] << 4) + 9; CN[(i + 2) * n + j] = (CN[(i + 2) * n + j] << 4) + 10;
                }
                if (j < n - 2)
                {
                    C[cnt++] = new constraint { x = i * n + j, y = i * n + j + 2, z = 4f / (n - 1) };
                    CN[i * n + j] = (CN[i * n + j] << 4) + 11; CN[i * n + j + 2] = (CN[i * n + j + 2] << 4) + 12;
                }
            }
    }
    
    void FixedUpdate()
    {
        // Get Vertices
        P = mesh.vertices; t = Time.realtimeSinceStartup;

        switch (type)
        {
            case Type.Semi_Implict:
                Semi_Implict();
                break;
            case Type.Semi_Implict_Burst:
                Semi_ImplictBurst(thread);
                break;
            case Type.Semi_Implict_Cpp:
                Semi_Implict_Cpp(P, P1, P2, transform.localToWorldMatrix, transform.worldToLocalMatrix, g, k, iter, dt, damping);
                break;
            case Type.Semi_Implict_CppThread:
                Semi_Implict_CppThread(P, P1, P2, transform.localToWorldMatrix, transform.worldToLocalMatrix, g, k, iter, thread, dt, damping);
                break;
        }

        // Update Mesh
        mesh.vertices = P;  mesh.RecalculateNormals();
        texts[0].text = (int)(1000 * (Time.realtimeSinceStartup - t)) + " ms";
    }

    [BurstCompile]
    struct TransformBurst: IJobParallelFor
    {
        [ReadOnly] public Matrix4x4 mat;
        public NativeArray<float3> P;

        public void Execute(int i)
        {
            P[i] = mat.MultiplyPoint3x4(P[i]);
        }
    }

    void Semi_Implict()
    {
        for (int i = 0; i < n * n; i++) P2[i] = transform.TransformPoint(P[i]);

        float t = dt / iter, t2 = t * t, d = Mathf.Pow(damping, 1f / iter);
        for (int i = 0; i < iter; i++) 
        {
            // Constrant
            for (int j = 0; j < C.Length; j++)
            {
                int x = C[j].x, y = C[j].y; float z = C[j].z;
                Vector3 f = k * ((P2[x] - P2[y]).magnitude - z) * (P2[x] - P2[y]).normalized;
                F[x] -= f; F[y] += f;
            }

            // Gravity & Update & Calculate Local Position
            for (int j = 0; j < n * n; j++)
            {
                if (j != n - 1 && j != n * n - 1) P2[j] += d * (P2[j] - P1[j]) + t2 * F[j];
                P1[j] = P2[j]; F[j] = g;
            }
        }

        for (int i = 0; i < n * n; i++) P[i] = transform.InverseTransformPoint(P2[i]);
    }

    void Semi_ImplictBurst(int thr)
    {
        MemCpy(P, P2N);
        handle = new TransformBurst() { mat = transform.localToWorldMatrix, P = P2N }.Schedule(n * n, n * n / thr); handle.Complete(); 

        float t = dt / iter, t2 = t * t, d = Mathf.Pow(damping, 1f / iter);
        for (int i = 0; i < iter; i++)
        {
            handle = new HandleConstaintBurst() { F = FN, P2 = P2N, C = CN, k = k }.Schedule(n * n, n * n / thr); handle.Complete();
            handle = new UpdatePositionBurst() { F = FN, P1 = P1N, P2 = P2N, g = g, d = d, t2 = t2 }.Schedule(n * n, n * n / thr); handle.Complete();
        }
        
        handle = new TransformBurst() { mat = transform.worldToLocalMatrix, P = P2N }.Schedule(n * n, n * n / thr); handle.Complete();
        MemCpy(P2N, P);
    }

    [BurstCompile]
    struct HandleConstaintBurst : IJobParallelFor
    {
        public NativeArray<float3> F;
        [ReadOnly] public NativeArray<float3> P2;
        [ReadOnly] public NativeArray<long> C;
        [ReadOnly] public float k;

        public void Execute(int i)
        {
            // Unpack C[i]
            long c = C[i]; float3 p;
            while (c > 0)
            {
                switch (c % 16)
                {
                    case 1: p = P2[i] - P2[i + n]; F[i] -= k * (math.length(p) - 2f / (n - 1)) * math.normalize(p); break;
                    case 2: p = P2[i] - P2[i - n]; F[i] -= k * (math.length(p) - 2f / (n - 1)) * math.normalize(p); break;
                    case 3: p = P2[i] - P2[i + 1]; F[i] -= k * (math.length(p) - 2f / (n - 1)) * math.normalize(p); break;
                    case 4: p = P2[i] - P2[i - 1]; F[i] -= k * (math.length(p) - 2f / (n - 1)) * math.normalize(p); break;
                    case 5: p = P2[i] - P2[i + n + 1]; F[i] -= k * (math.length(p) - 2.828f / (n - 1)) * math.normalize(p); break;
                    case 6: p = P2[i] - P2[i - n - 1]; F[i] -= k * (math.length(p) - 2.828f / (n - 1)) * math.normalize(p); break;
                    case 7: p = P2[i] - P2[i + n - 1]; F[i] -= k * (math.length(p) - 2.828f / (n - 1)) * math.normalize(p); break;
                    case 8: p = P2[i] - P2[i - n + 1]; F[i] -= k * (math.length(p) - 2.828f / (n - 1)) * math.normalize(p); break;
                    case 9: p = P2[i] - P2[i + 2 * n]; F[i] -= k * (math.length(p) - 4f / (n - 1)) * math.normalize(p); break;
                    case 10: p = P2[i] - P2[i - 2 * n]; F[i] -= k * (math.length(p) - 4f / (n - 1)) * math.normalize(p); break;
                    case 11: p = P2[i] - P2[i + 2]; F[i] -= k * (math.length(p) - 4f / (n - 1)) * math.normalize(p); break;
                    case 12: p = P2[i] - P2[i - 2]; F[i] -= k * (math.length(p) - 4f / (n - 1)) * math.normalize(p); break;
                }
                c >>= 4;
            }
        }
    }

    [BurstCompile]
    struct UpdatePositionBurst : IJobParallelFor
    {
        public NativeArray<float3> F, P1, P2;
        [ReadOnly] public float d, t2; [ReadOnly] public float3 g;

        public void Execute(int i)
        {
            if (i != n - 1 && i != n * n - 1) P2[i] += d * (P2[i] - P1[i]) + t2 * F[i];
            P1[i] = P2[i]; F[i] = g;
        }
    }

    unsafe static void MemCpy<SRC, DST>(SRC[] src, NativeArray<DST> dst) where SRC : struct where DST : struct
    {
        int srcSize = src.Length * UnsafeUtility.SizeOf<SRC>();
        int dstSize = dst.Length * UnsafeUtility.SizeOf<DST>();
        void* srcPtr = UnsafeUtility.PinGCArrayAndGetDataAddress(src, out ulong handle);
        void* dstPtr = NativeArrayUnsafeUtility.GetUnsafeReadOnlyPtr(dst);
        UnsafeUtility.MemCpy(destination: dstPtr, source: srcPtr, size: srcSize);
        UnsafeUtility.ReleaseGCObject(handle);
    }

    unsafe static void MemCpy<SRC, DST>(NativeArray<SRC> src, DST[] dst) where SRC : struct where DST : struct
    {
        int srcSize = src.Length * UnsafeUtility.SizeOf<SRC>();
        int dstSize = dst.Length * UnsafeUtility.SizeOf<DST>();
        void* srcPtr = NativeArrayUnsafeUtility.GetUnsafeReadOnlyPtr(src);
        void* dstPtr = UnsafeUtility.PinGCArrayAndGetDataAddress(dst, out ulong handle);
        UnsafeUtility.MemCpy(destination: dstPtr, source: srcPtr, size: srcSize);
        UnsafeUtility.ReleaseGCObject(handle);
    }

    public void ChangeSlide(int num)
    {
        switch (num)
        {
            case 0:
                texts[1].text = "k: " + (k = (int)sliders[0].value); break;
            case 1:
                texts[2].text = "dt: " + Math.Round(dt = sliders[1].value, 3); break;
            case 2:
                texts[3].text = "iter: " + (iter = (int)sliders[2].value); break;
        }
    }

    public void ChangeType(int num)
    {
        type = (Type)num;
        Start();
    }
}
