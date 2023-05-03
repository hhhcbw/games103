using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;

public class implicit_model : MonoBehaviour
{
    float t = 0.0333f;
    float mass = 1;
    float damping = 0.99f;
    float rho = 0.995f;
    float spring_k = 8000;
    Vector3 gravity = new Vector3(0, -9.8f, 0);
    int[] E;
    float[] L;
    Vector3[] V;

    // Start is called before the first frame update
    void Start()
    {
        Mesh mesh = GetComponent<MeshFilter>().mesh;

        //Resize the mesh.
        int n = 21;
        Vector3[] X = new Vector3[n * n];
        Vector2[] UV = new Vector2[n * n];
        int[] triangles = new int[(n - 1) * (n - 1) * 6];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++)
            {
                X[j * n + i] = new Vector3(5 - 10.0f * i / (n - 1), 0, 5 - 10.0f * j / (n - 1));
                UV[j * n + i] = new Vector3(i / (n - 1.0f), j / (n - 1.0f));
            }
        int t = 0;
        for (int j = 0; j < n - 1; j++)
            for (int i = 0; i < n - 1; i++)
            {
                triangles[t * 6 + 0] = j * n + i;
                triangles[t * 6 + 1] = j * n + i + 1;
                triangles[t * 6 + 2] = (j + 1) * n + i + 1;
                triangles[t * 6 + 3] = j * n + i;
                triangles[t * 6 + 4] = (j + 1) * n + i + 1;
                triangles[t * 6 + 5] = (j + 1) * n + i;
                t++;
            }
        mesh.vertices = X;
        mesh.triangles = triangles;
        mesh.uv = UV;
        mesh.RecalculateNormals();


        //Construct the original E
        int[] _E = new int[triangles.Length * 2];
        for (int i = 0; i < triangles.Length; i += 3)
        {
            _E[i * 2 + 0] = triangles[i + 0];
            _E[i * 2 + 1] = triangles[i + 1];
            _E[i * 2 + 2] = triangles[i + 1];
            _E[i * 2 + 3] = triangles[i + 2];
            _E[i * 2 + 4] = triangles[i + 2];
            _E[i * 2 + 5] = triangles[i + 0];
        }
        //Reorder the original edge list
        for (int i = 0; i < _E.Length; i += 2)
            if (_E[i] > _E[i + 1])
                Swap(ref _E[i], ref _E[i + 1]);
        //Sort the original edge list using quicksort
        Quick_Sort(ref _E, 0, _E.Length / 2 - 1);

        int e_number = 0;
        for (int i = 0; i < _E.Length; i += 2)
            if (i == 0 || _E[i + 0] != _E[i - 2] || _E[i + 1] != _E[i - 1])
                e_number++;

        E = new int[e_number * 2];
        for (int i = 0, e = 0; i < _E.Length; i += 2)
            if (i == 0 || _E[i + 0] != _E[i - 2] || _E[i + 1] != _E[i - 1])
            {
                E[e * 2 + 0] = _E[i + 0];
                E[e * 2 + 1] = _E[i + 1];
                e++;
            }

        L = new float[E.Length / 2];
        for (int e = 0; e < E.Length / 2; e++)
        {
            int v0 = E[e * 2 + 0];
            int v1 = E[e * 2 + 1];
            L[e] = (X[v0] - X[v1]).magnitude;
        }

        V = new Vector3[X.Length];
        for (int i = 0; i < V.Length; i++)
            V[i] = new Vector3(0, 0, 0);
    }

    void Quick_Sort(ref int[] a, int l, int r)
    {
        int j;
        if (l < r)
        {
            j = Quick_Sort_Partition(ref a, l, r);
            Quick_Sort(ref a, l, j - 1);
            Quick_Sort(ref a, j + 1, r);
        }
    }

    int Quick_Sort_Partition(ref int[] a, int l, int r)
    {
        int pivot_0, pivot_1, i, j;
        pivot_0 = a[l * 2 + 0];
        pivot_1 = a[l * 2 + 1];
        i = l;
        j = r + 1;
        while (true)
        {
            do ++i; while (i <= r && (a[i * 2] < pivot_0 || a[i * 2] == pivot_0 && a[i * 2 + 1] <= pivot_1));
            do --j; while (a[j * 2] > pivot_0 || a[j * 2] == pivot_0 && a[j * 2 + 1] > pivot_1);
            if (i >= j) break;
            Swap(ref a[i * 2], ref a[j * 2]);
            Swap(ref a[i * 2 + 1], ref a[j * 2 + 1]);
        }
        Swap(ref a[l * 2 + 0], ref a[j * 2 + 0]);
        Swap(ref a[l * 2 + 1], ref a[j * 2 + 1]);
        return j;
    }

    void Swap(ref int a, ref int b)
    {
        int temp = a;
        a = b;
        b = temp;
    }

    void Collision_Handling()
    {
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] X = mesh.vertices;

        //For every vertex, detect collision and apply impulse if needed.
        float radius = 2.7f;
        GameObject sphere = GameObject.Find("Sphere");
        Vector3 center = sphere.transform.position;
        for (int i = 0; i < X.Length; i++)
        {
            if (i == 0 || i == 20)
                continue;
            Vector3 d = X[i] - center;
            if (d.magnitude < radius)
            {
                Vector3 A = center + radius * d.normalized;
                V[i] = V[i] + (A - X[i]) / t;
                X[i] = A;
            }
        }

        mesh.vertices = X;
    }

    void Get_Gradient(Vector3[] X, Vector3[] X_hat, float t, Vector3[] G)
    {
        //Momentum and Gravity.
        for (int i = 0; i < X.Length; i++)
        {
            G[i] = mass * (X[i] - X_hat[i]) / (t * t) - mass * gravity;
        }
        //Spring Force.
        for (int e = 0; e < L.Length; e++)
        {
            int i = E[e * 2];
            int j = E[e * 2 + 1];
            G[i] = G[i] + spring_k * (1 - L[e] / (X[i] - X[j]).magnitude) * (X[i] - X[j]);
            G[j] = G[j] - spring_k * (1 - L[e] / (X[i] - X[j]).magnitude) * (X[i] - X[j]);
        }
    }

    // Update is called once per frame
    void Update()
    {
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] X = mesh.vertices;
        Vector3[] last_X = new Vector3[X.Length];
        Vector3[] X_hat = new Vector3[X.Length];
        Vector3[] G = new Vector3[X.Length];

        //Initial Setup.
        for (int i = 0; i < X.Length; i++)
        {
            if (i == 0 || i == 20) continue;
            V[i] = V[i] * damping;
            X_hat[i] = X[i] + t * V[i];
            X[i] = X_hat[i];
        }

        float omega = 1.0f;
        for (int k = 0; k < 32; k++)
        {
            //if (k == 0) omega = 1.0f;
            //else if (k == 1) omega = 2.0f / (2.0f - rho * rho);
            //else omega = 4.0f / (4.0f - rho * rho * omega);

            Get_Gradient(X, X_hat, t, G);

            //Update X by gradient.
            for (int i = 0; i < X.Length; i++)
            {
                if (i == 0 || i == 20) continue;
                Vector3 new_x = omega * (X[i] + (1.0f / (mass / (t * t) + 4.0f * spring_k)) * -G[i]) + (1.0f-omega) * last_X[i];
                last_X[i] = X[i];
                X[i] = new_x;
            }
        }

        //Finishing.
        for (int i = 0; i < X.Length; i++)
            V[i] += (X[i] - X_hat[i]) / t;
        mesh.vertices = X;

        Collision_Handling();
        mesh.RecalculateNormals();
    }
}
