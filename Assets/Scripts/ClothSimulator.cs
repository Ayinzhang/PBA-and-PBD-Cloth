using System;
using System.Runtime.InteropServices;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.Jobs;
using Unity.Jobs; // Define [IJob] [IJobParallelFor]
using Unity.Burst; // Define [BurstCompile]
using Unity.Collections; // Define [NativeArray]
using Unity.Mathematics; // SIMD Math Library
using Unity.Collections.LowLevel.Unsafe; // Quick Memory Copy

public class ClothSimulator : MonoBehaviour
{
    public enum Type { Semi_Implict, Semi_Implict_Burst, Semi_Implict_Cpp, Semi_Implict_CppThread, 
        Implict, Implict_Burst, PBD, PBD_Burst, XPBD, XPBD_Burst}
    public Type type = Type.Semi_Implict; public int k = 2000, iter = 32, thread = 64; public float dt = 1f / 30, damping = 0.95f;

    Mesh mesh; Text[] texts; Slider[] sliders; Vector3 g = new Vector3(0, -9.8f, 0); float t, rho = 0.995f; const int n = 21;
    Vector3[] P = new Vector3[n * n], P1 = new Vector3[n * n], P2 = new Vector3[n * n], F = new Vector3[n * n], G = new Vector3[n * n],
    V = new Vector3[n * n], SX = new Vector3[n * n]; Vector2[] UV = new Vector2[n * n]; int[] T = new int[(n - 1) * (n - 1) * 6]; 
    float[] W = new float[n * n]; int[] SN = new int[n * n];
    struct constraint { public int x, y; public float z; } constraint[] C = new constraint[5 * n * n - 8 * n + 1];

    JobHandle handle; NativeArray<float3> FN = new NativeArray<float3>(n * n, Allocator.Persistent), 
    GN = new NativeArray<float3>(n * n, Allocator.Persistent), PN = new NativeArray<float3>(n * n, Allocator.Persistent), 
    P1N = new NativeArray<float3>(n * n, Allocator.Persistent), P2N = new NativeArray<float3>(n * n, Allocator.Persistent), 
    VN = new NativeArray<float3>(n * n, Allocator.Persistent), SXN = new NativeArray<float3>(n * n, Allocator.Persistent);
    NativeArray<long> CN = new NativeArray<long>(n * n, Allocator.Persistent);
    NativeArray<float> WN = new NativeArray<float>(n * n, Allocator.Persistent);
    NativeArray<int> SNN = new NativeArray<int>(n * n, Allocator.Persistent);

    [DllImport("ClothSimulator", CallingConvention = CallingConvention.Cdecl)]
    static extern void Semi_Implict_Cpp( [In, Out] Vector3[] P, [In, Out] Vector3[] P1, [In, Out] Vector3[] P2, [In, Out] Vector3[] V, 
        Matrix4x4 mat, Matrix4x4 imat, Vector3 g, int k, int iter, float dt, float damping);

    [DllImport("ClothSimulator", CallingConvention = CallingConvention.Cdecl)]
    static extern void Semi_Implict_CppThread([In, Out] Vector3[] P, [In, Out] Vector3[] P1, [In, Out] Vector3[] P2, [In, Out] Vector3[] V,
        Matrix4x4 mat, Matrix4x4 imat, Vector3 g, int k, int iter, int thr, float dt, float damping);

    void Start()
    {
        transform.position = Vector3.zero;
        mesh = GetComponent<MeshFilter>().mesh;
        texts = GameObject.Find("Canvas").GetComponentsInChildren<Text>();
        sliders = GameObject.Find("Canvas").GetComponentsInChildren<Slider>();

        // Resize the mesh (Blender type)
        for (int i = 0; i < n * n; i++) if (i != n - 1 && i != n * n - 1) W[i] = WN[i] = 1;

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                P[i * n + j] = PN[i * n + j] = new Vector3(2f * i / (n - 1) - 1, 2f * j / (n - 1) - 1, 0);
                P1[i * n + j] = P1N[i * n + j] = P2[i * n + j] = P2N[i * n + j] = transform.TransformPoint(P[i * n + j]);
                F[i * n + j] = FN[i * n + j] = G[i * n + j] = GN[i * n + j] = Vector3.zero;
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
        cnt = 0;
        for (int i = 0; i < n * n; i++) CN[i] = 0;
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
                Semi_Implict_Cpp(P, P1, P2, V, transform.localToWorldMatrix, transform.worldToLocalMatrix, g, k, iter, dt, damping);
                break;
            case Type.Semi_Implict_CppThread:
                Semi_Implict_CppThread(P, P1, P2, V, transform.localToWorldMatrix, transform.worldToLocalMatrix, g, k, iter, thread, dt, damping);
                break;
            case Type.Implict:
                Implict();
                break;
            case Type.Implict_Burst:
                ImplictBurst(thread);
                break;
            case Type.PBD:
                PBD();
                break;
            case Type.PBD_Burst:
                PBDBurst(thread);
                break;
            case Type.XPBD:
                XPBD();
                break;
            case Type.XPBD_Burst:
                XPBDBurst(thread);
                break;
        }

        // Update Mesh
        mesh.vertices = P;  mesh.RecalculateNormals();
        texts[0].text = String.Format("{0:N2} ms", 1000 * (Time.realtimeSinceStartup - t));
    }

    void Semi_Implict()
    {
        for (int i = 0; i < n * n; i++) 
        { 
            P2[i] = transform.TransformPoint(P[i]); 
            V[i] += (P2[i] - P1[i]) / dt;
        }

        float idt = dt / iter, d = Mathf.Pow(damping, 1f / iter);
        for (int i = 0; i < iter; i++)
        {
            // Constrant
            for (int j = 0; j < 5 * n * n - 8 * n + 1; j++)
            {
                int x = C[j].x, y = C[j].y; float z = C[j].z;
                Vector3 f = k * (1 - z / (P2[x] - P2[y]).magnitude) * (P2[x] - P2[y]);
                F[x] -= f; F[y] += f;
            }

            // Gravity & Update & Calculate Local Position
            for (int j = 0; j < n * n; j++)
            {
                V[j] = d * V[j] + idt * F[j];
                if (W[j] != 0) P2[j] += idt * V[j];
                P1[j] = P2[j]; F[j] = g;
            }
        }

        for (int i = 0; i < n * n; i++) P[i] = transform.InverseTransformPoint(P2[i]);
    }

    void Semi_ImplictBurst(int thr)
    {
        MemCpy(P, P2N);
        handle = new SemiImplictTransformInBurst() { mat = transform.localToWorldMatrix, P = P2N, P1 = P1N, V = VN, dt = dt }.Schedule(n * n, n * n / thr); handle.Complete();

        float idt = dt / iter, d = Mathf.Pow(damping, 1f / iter);
        for (int i = 0; i < iter; i++)
        {
            handle = new HandleConstaintBurst() { F = FN, P2 = P2N, C = CN, k = k }.Schedule(n * n, n * n / thr); handle.Complete();
            handle = new UpdatePositionBurst() { V = VN, F = FN, W = WN, P1 = P1N, P2 = P2N, g = g, d = d, idt = idt }.Schedule(n * n, n * n / thr); handle.Complete();
        }

        handle = new TransformBurst() { mat = transform.worldToLocalMatrix, P = P2N }.Schedule(n * n, n * n / thr); handle.Complete();
        MemCpy(P2N, P);
    }

    void Implict()
    {
        for (int i = 0; i < n * n; i++)
        {
            V[i] += ((P2[i] = transform.TransformPoint(P[i])) - P1[i]) / dt;
            if (W[i] != 0) P[i] = P2[i] += (dt * (V[i] *= damping));
        }

        float invdt2 = 1 / dt / dt, factor = 1 / (invdt2 + 4 * k), rho2 = rho * rho, w = 1;
        for (int i = 0; i < iter; i++)
        {
            //Momentum and Gravity.
            for (int j = 0; j < n * n; j++) G[j] = invdt2 * (P2[j] - P[j]) - g;

            //Spring Force.
            for (int j = 0; j < 5 * n * n - 8 * n + 1; j++)
            {
                int x = C[j].x, y = C[j].y; float z = C[j].z;
                Vector3 grad = k * (1 - z / (P2[x] - P2[y]).magnitude) * (P2[x] - P2[y]);
                G[x] += grad; G[y] -= grad;
            }

            if (i == 0) w = 1;
            else if (i == 1) w = 2 / (2 - rho2);
            else w = 4 / (4 - w * rho2);

            //Update P2 by gradient.
            for (int j = 0; j < n * n; j++)
            {
                if (W[j] == 0) continue;
                Vector3 old = P2[j];
                P2[j] = w * (P2[j] - factor * G[j]) + (1 - w) * P1[j];
                P1[j] = old;
            }
        }

        for (int i = 0; i < n * n; i++) 
        { 
            V[i] += (P2[i] - P[i]) / dt; 
            P[i] = transform.InverseTransformPoint(P2[i]); 
        }
    }

    void ImplictBurst(int thr)
    {
        MemCpy(P, P2N);
        handle = new ImplictTransformInBurst() { P = PN, P1 = P1N, P2 = P2N, W = WN, V = VN, 
            mat = transform.localToWorldMatrix, dt = dt, damping = damping }.Schedule(n * n, n * n / thr); handle.Complete();

        float invdt2 = 1 / dt / dt, factor = 1 / (invdt2 + 4 * k);
        for (int i = 0; i < iter; i++)
        { 
            handle = new ImplictCalculateGradientBurst(){ P = PN, P2 = P2N, C = CN, G = GN, W = WN,
                g = g, k = k, invdt2 = invdt2 }.Schedule(n * n, n * n / thr); handle.Complete();
            handle = new ImplictUpdatePositionBurst() { P1 = P1N, P2 = P2N, G = GN, W = WN,
                factor = factor }.Schedule(n * n, n * n / thr); handle.Complete();
        }

        handle = new ImplictTransformInBurst() { P = PN, P1 = P1N, P2 = P2N, V = VN, W = WN,
            mat = transform.worldToLocalMatrix, dt = dt }.Schedule(n * n, n * n / thr); handle.Complete();
        MemCpy(P2N, P);
    }

    void PBD()
    {
        for (int i = 0; i < n * n; i++) P2[i] = transform.TransformPoint(P[i]);

        // Free Update
        for (int i = 0; i < n * n; i++)
        {
            if (W[i] == 0) continue;
            V[i] = damping * V[i] + dt * g;
            P2[i] += dt * V[i] + damping * (P2[i] - P1[i]);
        }

        for (int i = 0; i < iter; i++)
        {
            // Strain_Limit
            for (int j = 0; j < 5 * n * n - 8 * n + 1; j++)
            {
                int x = C[j].x, y = C[j].y; float z = C[j].z;
                Vector3 p = z * (P2[x] - P2[y]).normalized;
                SX[x] += 0.5f * (P2[x] + P2[y] + p); SN[x]++; SX[y] += 0.5f * (P2[x] + P2[y] - p); SN[y]++;
            }
            for (int j = 0; j < n * n; j++)
            {
                if (W[j] == 0) continue;
                P1[j] = P2[j];
                P2[j] = (0.2f * P2[j] + SX[j]) / (0.2f + SN[j]);
                V[j] += (P2[j] - P1[j]) / dt;
                SX[j] = Vector3.zero; SN[j] = 0;
            }
        }
        
        for (int i = 0; i < n * n; i++) P[i] = transform.InverseTransformPoint(P1[i] = P2[i]);
    }

    void PBDBurst(int thr)
    {
        MemCpy(P, P2N);
        handle = new PBDTransformInBurst() { P1 = P1N, P2 = P2N, W = WN, V = VN,
            mat = transform.localToWorldMatrix, g = g, dt = dt, damping = damping }.Schedule(n * n, n * n / thr); handle.Complete();
        for (int i = 0; i < iter; i++)
        {
            handle = new PBDCalculateGradientBurst() { P2 = P2N, SX = SXN, SN = SNN, C = CN
                }.Schedule(n * n, n * n / thr); handle.Complete();
            handle = new PBDUpdatePositionBurst() { P1 = P1N, P2 = P2N, SX = SXN, SN = SNN, W = WN, V = VN,
                dt = dt }.Schedule(n * n, n * n / thr); handle.Complete();
        }
        handle = new PBDTransformOutBurst() { P1 = P1N, P2 = P2N,
            mat = transform.worldToLocalMatrix }.Schedule(n * n, n * n / thr); handle.Complete();
        MemCpy(P2N, P);
    }

    void XPBD()
    {
        float idt = dt / iter, a = 1f / k / idt / idt, d = Mathf.Pow(damping, 1f / iter);
        for (int i = 0; i < n * n; i++) { V[i] += ((P2[i] = transform.TransformPoint(P[i])) - P1[i]) / idt; }

        for (int i = 0; i < iter; i++)
        {
            // Free Update
            for (int j = 0; j < n * n; j++)
            {
                if (W[j] == 0) continue;
                V[j] = d * V[j] + idt * g;
                P2[j] += idt * V[j];
            }

            // Solve Constraint (Gauss Seidel)
            for (int j = 0; j < 5 * n * n - 8 * n + 1; j++)
            {
                int x = C[j].x, y = C[j].y; float z = C[j].z, w = W[x] + W[y];
                Vector3 grad = -(1 - z / (P2[x] - P2[y]).magnitude) * (P2[x] - P2[y]) / (a + w);
                P2[x] += grad * W[x]; P2[y] -= grad * W[y];
            }

            // Record Velocity
            for (int j = 0; j < n * n; j++) { V[j] = (P2[j] - P1[j]) / idt; P1[j] = P2[j]; }
        }

        for (int i = 0; i < n * n; i++) P[i] = transform.InverseTransformPoint(P2[i]);
    }

    void XPBDBurst(int thr)
    {
        MemCpy(P, P2N);

        float idt = dt / iter, a = 1f / k / idt / idt, d = Mathf.Pow(damping, 1f / iter);
        handle = new XPBDTransformBurst() { P1 = P1N, P2 = P2N, V = VN,
        mat = transform.localToWorldMatrix, idt = idt }.Schedule(n * n, n * n / thr); handle.Complete();

        for (int i = 0; i < iter; i++)
        {
            handle = new XPBDFreeUpdateBurst() { P2 = P2N, V = VN, W = WN,
                g = g, d = damping, idt = idt }.Schedule(n * n, n * n / thr); handle.Complete();
            handle = new XPBDCalculateGradientBurst() { P2 = P2N, G = GN, C = CN, W = WN,
                g = g, a = a, d = damping, idt = idt }.Schedule(n * n, n * n / thr); handle.Complete();
            handle = new XPBDUpdatePositionBurst() { P1 = P1N, P2 = P2N, V = VN, G = GN,
                idt = idt }.Schedule(n * n, n * n / thr); handle.Complete();
        }

        handle = new TransformBurst() { P = P2N, mat = transform.worldToLocalMatrix }.Schedule(n * n, n * n / thr); handle.Complete();
        MemCpy(P2N, P);
    }

    [BurstCompile]
    struct TransformBurst: IJobParallelFor
    {
        public NativeArray<float3> P;
        [ReadOnly] public Matrix4x4 mat;

        public void Execute(int i)
        {
            P[i] = mat.MultiplyPoint3x4(P[i]);
        }
    }

    [BurstCompile]
    struct SemiImplictTransformInBurst : IJobParallelFor
    {
        public NativeArray<float3> P, P1, V;
        [ReadOnly] public Matrix4x4 mat;
        [ReadOnly] public float dt;

        public void Execute(int i)
        {
            P[i] = mat.MultiplyPoint3x4(P[i]);
            V[i] = (P[i] - P1[i]) / dt;
        }
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
        public NativeArray<float3> V, F, P1, P2;
        [ReadOnly] public NativeArray<float> W;
        [ReadOnly] public float d, idt; [ReadOnly] public float3 g;

        public void Execute(int i)
        {
            V[i] = d * V[i] + idt * F[i];
            if (W[i] != 0) P2[i] += idt * V[i];
            P1[i] = P2[i]; F[i] = g;
        }
    }

    [BurstCompile]
    struct ImplictTransformInBurst : IJobParallelFor
    {
        public NativeArray<float3> P, P1, P2, V;
        [ReadOnly] public NativeArray<float> W;
        [ReadOnly] public Matrix4x4 mat;
        [ReadOnly] public float dt, damping;

        public void Execute(int i)
        {
            P2[i] = mat.MultiplyPoint3x4(P2[i]);
            V[i] += (P2[i] - P1[i]) / dt;
            if (W[i] != 0) P[i] = P2[i] += (dt * (V[i] *= damping));
        }
    }

    [BurstCompile]
    struct ImplictCalculateGradientBurst : IJobParallelFor
    {
        public NativeArray<float3> G;
        [ReadOnly] public NativeArray<float3> P, P2;
        [ReadOnly] public NativeArray<long> C;
        [ReadOnly] public NativeArray<float> W;
        [ReadOnly] public float3 g;
        [ReadOnly] public float k, invdt2;

        public void Execute(int i)
        {
            G[i] = invdt2 * (P2[i] - P[i]) - g;
            long c = C[i]; float3 p;
            while (c > 0)
            {
                switch (c % 16)
                {
                    case 1: p = P2[i] - P2[i + n]; G[i] += k * (1 - 2f / (n - 1) / math.length(p)) * p; break;
                    case 2: p = P2[i] - P2[i - n]; G[i] += k * (1 - 2f / (n - 1) / math.length(p)) * p; break;
                    case 3: p = P2[i] - P2[i + 1]; G[i] += k * (1 - 2f / (n - 1) / math.length(p)) * p; break;
                    case 4: p = P2[i] - P2[i - 1]; G[i] += k * (1 - 2f / (n - 1) / math.length(p)) * p; break;
                    case 5: p = P2[i] - P2[i + n + 1]; G[i] += k * (1 - 2.828f / (n - 1) / math.length(p)) * p; break;
                    case 6: p = P2[i] - P2[i - n - 1]; G[i] += k * (1 - 2.828f / (n - 1) / math.length(p)) * p; break;
                    case 7: p = P2[i] - P2[i + n - 1]; G[i] += k * (1 - 2.828f / (n - 1) / math.length(p)) * p; break;
                    case 8: p = P2[i] - P2[i - n + 1]; G[i] += k * (1 - 2.828f / (n - 1) / math.length(p)) * p; break;
                    case 9: p = P2[i] - P2[i + 2 * n]; G[i] += k * (1 - 4f / (n - 1) / math.length(p)) * p; break;
                    case 10: p = P2[i] - P2[i - 2 * n]; G[i] += k * (1 - 4f / (n - 1) / math.length(p)) * p; break;
                    case 11: p = P2[i] - P2[i + 2]; G[i] += k * (1 - 4f / (n - 1) / math.length(p)) * p; break;
                    case 12: p = P2[i] - P2[i - 2]; G[i] += k * (1 - 4f / (n - 1) / math.length(p)) * p; break;
                }
                c >>= 4;
            }
        }
    }

    [BurstCompile]
    struct ImplictUpdatePositionBurst : IJobParallelFor
    {
        public NativeArray<float3> P1, P2, G;
        [ReadOnly] public NativeArray<float> W;
        [ReadOnly] public float factor;

        public void Execute(int i)
        {
            P1[i] = P2[i];
            if (W[i] != 0) P2[i] -= factor * G[i];
        }
    }

    [BurstCompile]
    struct ImplictTransformOutBurst : IJobParallelFor
    {
        public NativeArray<float3> P2, V;
        [ReadOnly] NativeArray<float3> P;
        [ReadOnly] public Matrix4x4 mat;
        [ReadOnly] public float dt;

        public void Execute(int i)
        {
            V[i] += (P2[i] - P[i]) / dt;
            P2[i] = mat.MultiplyPoint3x4(P2[i]);
        }
    }

    [BurstCompile]
    struct PBDTransformInBurst : IJobParallelFor
    {
        public NativeArray<float3> P2, V;
        [ReadOnly] public NativeArray<float3> P1;
        [ReadOnly] public NativeArray<float> W;
        [ReadOnly] public Matrix4x4 mat;
        [ReadOnly] public float3 g;
        [ReadOnly] public float dt, damping;

        public void Execute(int i)
        {
            P2[i] = mat.MultiplyPoint3x4(P2[i]);
            if (W[i] == 0) return;
            V[i] = damping * V[i] + dt * g;
            P2[i] += dt * V[i] + damping * (P2[i] - P1[i]);
        }
    }

    [BurstCompile]
    struct PBDCalculateGradientBurst : IJobParallelFor
    {
        public NativeArray<float3> SX;
        public NativeArray<int> SN;
        [ReadOnly] public NativeArray<float3> P2;
        [ReadOnly] public NativeArray<long> C;

        public void Execute(int i)
        {
            long c = C[i];
            while (c > 0)
            {
                switch (c % 16)
                {
                    case 1: SX[i] += 0.5f * (P2[i] + P2[i + n] + 2f / (n - 1) * math.normalize(P2[i] - P2[i + n])); break;
                    case 2: SX[i] += 0.5f * (P2[i] + P2[i - n] + 2f / (n - 1) * math.normalize(P2[i] - P2[i - n])); break;
                    case 3: SX[i] += 0.5f * (P2[i] + P2[i + 1] + 2f / (n - 1) * math.normalize(P2[i] - P2[i + 1])); break;
                    case 4: SX[i] += 0.5f * (P2[i] + P2[i - 1] + 2f / (n - 1) * math.normalize(P2[i] - P2[i - 1])); break;
                    case 5: SX[i] += 0.5f * (P2[i] + P2[i + n + 1] + 2.828f / (n - 1) * math.normalize(P2[i] - P2[i + n + 1])); break;
                    case 6: SX[i] += 0.5f * (P2[i] + P2[i - n - 1] + 2.828f / (n - 1) * math.normalize(P2[i] - P2[i - n - 1])); break;
                    case 7: SX[i] += 0.5f * (P2[i] + P2[i + n - 1] + 2.828f / (n - 1) * math.normalize(P2[i] - P2[i + n - 1])); break;
                    case 8: SX[i] += 0.5f * (P2[i] + P2[i - n + 1] + 2.828f / (n - 1) * math.normalize(P2[i] - P2[i - n + 1])); break;
                    case 9: SX[i] += 0.5f * (P2[i] + P2[i + 2 * n] + 4f / (n - 1) * math.normalize(P2[i] - P2[i + 2 * n])); break;
                    case 10: SX[i] += 0.5f * (P2[i] + P2[i - 2 * n] + 4f / (n - 1) * math.normalize(P2[i] - P2[i - 2 * n])); break;
                    case 11: SX[i] += 0.5f * (P2[i] + P2[i + 2] + 4f / (n - 1) * math.normalize(P2[i] - P2[i + 2])); break;
                    case 12: SX[i] += 0.5f * (P2[i] + P2[i - 2] + 4f / (n - 1) * math.normalize(P2[i] - P2[i - 2])); break;
                }
                c >>= 4; SN[i]++;
            }
        }
    }

    [BurstCompile]
    struct PBDUpdatePositionBurst : IJobParallelFor
    {
        public NativeArray<float3> P1, P2, V, SX;
        public NativeArray<int> SN;
        [ReadOnly] public NativeArray<float> W;
        [ReadOnly] public float dt;

        public void Execute(int i)
        {
            if (W[i] == 0) return;
            P1[i] = P2[i];
            P2[i] = (0.2f * P2[i] + SX[i]) / (0.2f + SN[i]);
            V[i] += (P2[i] - P1[i]) / dt;
            SX[i] = float3.zero; SN[i] = 0;
        }
    }

    [BurstCompile]
    struct PBDTransformOutBurst : IJobParallelFor
    {
        public NativeArray<float3> P1, P2;
        [ReadOnly] public Matrix4x4 mat;

        public void Execute(int i)
        {
            P1[i] = P2[i];
            P2[i] = mat.MultiplyPoint3x4(P2[i]);
        }
    }

    [BurstCompile]
    struct XPBDTransformBurst : IJobParallelFor
    {
        public NativeArray<float3> P2, V;
        [ReadOnly] public NativeArray<float3> P1;
        [ReadOnly] public Matrix4x4 mat;
        [ReadOnly] public float idt;

        public void Execute(int i)
        {
            P2[i] = mat.MultiplyPoint3x4(P2[i]);
            V[i] += (P2[i] - P1[i]) / idt;
        }
    }

    [BurstCompile]
    struct XPBDFreeUpdateBurst : IJobParallelFor
    {
        public NativeArray<float3> P2, V;
        [ReadOnly] public NativeArray<float> W;
        [ReadOnly] public float3 g;
        [ReadOnly] public float d, idt;

        public void Execute(int i)
        {
            if (W[i] != 0)
            {
                V[i] = d * V[i] + idt * g;
                P2[i] += idt * V[i];
            }
        }
    }

    [BurstCompile]
    struct XPBDCalculateGradientBurst : IJobParallelFor
    {
        public NativeArray<float3> G;
        [ReadOnly] public NativeArray<float3> P2;
        [ReadOnly] public NativeArray<long> C;
        [ReadOnly] public NativeArray<float> W;
        [ReadOnly] public float3 g;
        [ReadOnly] public float a, d, idt;

        public void Execute(int i)
        {
            float3 p; long c = C[i];
            while (c > 0)
            {
                switch (c % 16)
                {
                    case 1: p = P2[i] - P2[i + n]; G[i] -= p * W[i] * (1 - 2f / (n - 1) / math.length(p)) / (a + W[i] + W[i + n]); break;
                    case 2: p = P2[i] - P2[i - n]; G[i] -= p * W[i] * (1 - 2f / (n - 1) / math.length(p)) / (a + W[i] + W[i - n]); break;
                    case 3: p = P2[i] - P2[i + 1]; G[i] -= p * W[i] * (1 - 2f / (n - 1) / math.length(p)) / (a + W[i] + W[i + 1]); break;
                    case 4: p = P2[i] - P2[i - 1]; G[i] -= p * W[i] * (1 - 2f / (n - 1) / math.length(p)) / (a + W[i] + W[i - 1]); break;
                    case 5: p = P2[i] - P2[i + n + 1]; G[i] -= p * W[i] * (1 - 2.828f / (n - 1) / math.length(p)) / (a + W[i] + W[i + n + 1]); break;
                    case 6: p = P2[i] - P2[i - n - 1]; G[i] -= p * W[i] * (1 - 2.828f / (n - 1) / math.length(p)) / (a + W[i] + W[i - n - 1]); break;
                    case 7: p = P2[i] - P2[i + n - 1]; G[i] -= p * W[i] * (1 - 2.828f / (n - 1) / math.length(p)) / (a + W[i] + W[i + n - 1]); break;
                    case 8: p = P2[i] - P2[i - n + 1]; G[i] -= p * W[i] * (1 - 2.828f / (n - 1) / math.length(p)) / (a + W[i] + W[i - n + 1]); break;
                    case 9: p = P2[i] - P2[i + 2 * n]; G[i] -= p * W[i] * (1 - 4f / (n - 1) / math.length(p)) / (a + W[i] + W[i + 2 * n]); break;
                    case 10: p = P2[i] - P2[i - 2 * n]; G[i] -= p * W[i] * (1 - 4f / (n - 1) / math.length(p)) / (a + W[i] + W[i - 2 * n]); break;
                    case 11: p = P2[i] - P2[i + 2]; G[i] -= p * W[i] * (1 - 4f / (n - 1) / math.length(p)) / (a + W[i] + W[i + 2]); break;
                    case 12: p = P2[i] - P2[i - 2]; G[i] -= p * W[i] * (1 - 4f / (n - 1) / math.length(p)) / (a + W[i] + W[i - 2]); break;
                }
                c >>= 4;
            }
        }
    }

    [BurstCompile]
    struct XPBDUpdatePositionBurst : IJobParallelFor
    {
        public NativeArray<float3> P1, P2, V, G;
        [ReadOnly] public float idt;

        public void Execute(int i)
        {
            P2[i] += G[i]; 
            G[i] = float3.zero;
            V[i] = (P2[i] - P1[i]) / idt; 
            P1[i] = P2[i];
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