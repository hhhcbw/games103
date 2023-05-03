using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.IO;
using System.Xml.XPath;
using UnityEditor;

public class FVM : MonoBehaviour
{
	float dt 			= 0.003f;
    float mass 			= 1;
	float stiffness_0	= 20000.0f;
    float stiffness_1 	= 5000.0f;
    float damp			= 0.999f;
    float restitution = 0.5f;
    float friction = 0.5f;


    int[] 		Tet;
	int tet_number;			//The number of tetrahedra

    Vector3 gravity = new Vector3(0.0f, -9.8f, 0.0f);
    Vector3 P = new Vector3(0.0f, -3.0f, 0.0f);
    Vector3 N = new Vector3(0.0f, 1.0f, 0.0f);

	Vector3[] 	Force;
	Vector3[] 	V;
	Vector3[] 	X;

	int number;				//The number of vertices

	Matrix4x4[] inv_Dm;

	//For Laplacian smoothing.
	Vector3[]   V_sum;
	int[]		V_num;

	SVD svd = new SVD();

    // Start is called before the first frame update
    void Start()
    {
    	// FILO IO: Read the house model from files.
    	// The model is from Jonathan Schewchuk's Stellar lib.
    	{
    		string fileContent = File.ReadAllText("Assets/house2.ele");
    		string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);
    		
    		tet_number=int.Parse(Strings[0]);
        	Tet = new int[tet_number*4];

            // index of vertices in tetrahedral
    		for(int tet=0; tet<tet_number; tet++)
    		{
				Tet[tet*4+0]=int.Parse(Strings[tet*5+4])-1;
				Tet[tet*4+1]=int.Parse(Strings[tet*5+5])-1;
				Tet[tet*4+2]=int.Parse(Strings[tet*5+6])-1;
				Tet[tet*4+3]=int.Parse(Strings[tet*5+7])-1;
			}
    	}
    	{
			string fileContent = File.ReadAllText("Assets/house2.node");
    		string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);
    		number = int.Parse(Strings[0]);
    		X = new Vector3[number];
       		for(int i=0; i<number; i++)
       		{
       			X[i].x=float.Parse(Strings[i*5+5])*0.4f;
       			X[i].y=float.Parse(Strings[i*5+6])*0.4f;
       			X[i].z=float.Parse(Strings[i*5+7])*0.4f;
       		}
    		//Centralize the model.
	    	Vector3 center=Vector3.zero;
	    	for(int i=0; i<number; i++)		center+=X[i];
	    	center=center/number;
	    	for(int i=0; i<number; i++)
	    	{
	    		X[i]-=center;
                // swap y and z
	    		float temp=X[i].y;
	    		X[i].y=X[i].z;
	    		X[i].z=temp;
	    	}
		}
        /*tet_number=1;
        Tet = new int[tet_number*4];
        Tet[0]=0;
        Tet[1]=1;
        Tet[2]=2;
        Tet[3]=3;

        number=4;
        X = new Vector3[number];
        V = new Vector3[number];
        Force = new Vector3[number];
        X[0]= new Vector3(0, 0, 0);
        X[1]= new Vector3(1, 0, 0);
        X[2]= new Vector3(0, 1, 0);
        X[3]= new Vector3(0, 0, 1);*/


        //Create triangle mesh.
       	Vector3[] vertices = new Vector3[tet_number*12];
        int vertex_number=0;
        // one tetrahedral mesh has four triangle mesh
        for(int tet=0; tet<tet_number; tet++)
        {
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];

        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];

        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];

        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        }

        int[] triangles = new int[tet_number*12];
        for(int t=0; t<tet_number*4; t++)
        {
        	triangles[t*3+0]=t*3+0;
        	triangles[t*3+1]=t*3+1;
        	triangles[t*3+2]=t*3+2;
        }
        Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.vertices  = vertices;
		mesh.triangles = triangles;
		mesh.RecalculateNormals ();


		V 	  = new Vector3[number];
        Force = new Vector3[number];
        V_sum = new Vector3[number];
        V_num = new int[number];

        //TODO: Need to allocate and assign inv_Dm
        inv_Dm = new Matrix4x4[tet_number];
        for (int tet = 0; tet < tet_number; tet++) 
        {
            inv_Dm[tet] = Build_Edge_Matrix(tet).inverse;
        }
    }

    Matrix4x4 Build_Edge_Matrix(int tet)
    {
    	Matrix4x4 ret=Matrix4x4.zero;
    	//TODO: Need to build edge matrix here.

        for (int i=0; i<3; i++) 
        {
            ret[0, i] = X[Tet[tet * 4 + i + 1]].x - X[Tet[tet * 4 + 0]].x;
            ret[1, i] = X[Tet[tet * 4 + i + 1]].y - X[Tet[tet * 4 + 0]].y;
            ret[2, i] = X[Tet[tet * 4 + i + 1]].z - X[Tet[tet * 4 + 0]].z;
        }
        ret[3, 3] = 1.0f;

		return ret;
    }

    Matrix4x4 Matrix4x4_minus_Matrix4x4(Matrix4x4 m1, Matrix4x4 m2)
    {
        Matrix4x4 res = Matrix4x4.zero;
        for (int i=0; i<4;  i++)
            for (int j=0; j<4; j++)
            {
                res[i, j] = m1[i, j] - m2[i, j];
            }
        return res;
    }

    Matrix4x4 Matrix4x4_add_Matrix4x4(Matrix4x4 m1, Matrix4x4 m2)
    {
        Matrix4x4 res = Matrix4x4.zero;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                res[i, j] = m1[i, j] + m2[i, j];
            }
        return res;
    }

    Matrix4x4 Matrix4x4_mul_float(Matrix4x4 m, float x)
    {
        Matrix4x4 res = Matrix4x4.zero;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
            {
                res[i, j] = m[i, j] * x;
            }
        return res;
    }

    void Smooth_V(float w)
    {
        for (int i = 0; i < number; i++)
        {
            V_sum[i] = Vector3.zero;
            V_num[i] = 0;
        }

        for (int tet = 0; tet < tet_number; tet++)
        {
            Vector3 sum = V[Tet[tet * 4 + 0]] + V[Tet[tet * 4 + 1]] + V[Tet[tet * 4 + 2]] + V[Tet[tet * 4 + 3]];
            V_sum[Tet[tet * 4 + 0]] += sum - V[Tet[tet * 4 + 0]];
            V_sum[Tet[tet * 4 + 1]] += sum - V[Tet[tet * 4 + 1]];
            V_sum[Tet[tet * 4 + 2]] += sum - V[Tet[tet * 4 + 2]];
            V_sum[Tet[tet * 4 + 3]] += sum - V[Tet[tet * 4 + 3]];
            V_num[Tet[tet * 4 + 0]] += 3;
            V_num[Tet[tet * 4 + 1]] += 3;
            V_num[Tet[tet * 4 + 2]] += 3;
            V_num[Tet[tet * 4 + 3]] += 3;
        }

        for (int i = 0; i < number; i++)
        {
            V[i] = w * V[i] + (1.0f - w) * V_sum[i] / V_num[i];
        }
    }

    void _Update()
    {
    	// Jump up.
		if(Input.GetKeyDown(KeyCode.Space))
    	{
    		for(int i=0; i<number; i++)
    			V[i].y+=0.2f;
    	}

    	for(int i=0 ;i<number; i++)
    	{
            //TODO: Add gravity to Force.
            Force[i] = gravity * mass;
    	}

    	for(int tet=0; tet<tet_number; tet++)
    	{
            //TODO: Deformation Gradient
            Matrix4x4 F = Build_Edge_Matrix(tet) * inv_Dm[tet];
            //TODO: Green Strain
            Matrix4x4 G = Matrix4x4_mul_float(Matrix4x4_minus_Matrix4x4(F.transpose * F, Matrix4x4.identity), 0.5f);
            //TODO: Second PK Stress
            Matrix4x4 S = Matrix4x4.zero;
            float trace = G[0, 0] + G[1, 1] + G[2, 2];
            S = Matrix4x4_add_Matrix4x4(Matrix4x4_mul_float(G, 2.0f * stiffness_1), Matrix4x4_mul_float(Matrix4x4.identity, stiffness_0 * trace));
            //TODO: Elastic Force
            float volume = 1.0f / (inv_Dm[tet].determinant * 6);
            Matrix4x4 Elastic_force = Matrix4x4_mul_float(F * S * inv_Dm[tet].transpose, -volume);
            for (int k=0; k<3; k++)
            {
                Force[Tet[tet * 4 + k + 1]].x += Elastic_force[0, k];
                Force[Tet[tet * 4 + k + 1]].y += Elastic_force[1, k];
                Force[Tet[tet * 4 + k + 1]].z += Elastic_force[2, k];
            }

            // Elastic_force0 = -Elastic_force1-Elastic_force2-Elastic_force3
            Force[Tet[tet * 4 + 0]].x -= Elastic_force[0, 0] + Elastic_force[0, 1] + Elastic_force[0, 2];
            Force[Tet[tet * 4 + 0]].y -= Elastic_force[1, 0] + Elastic_force[1, 1] + Elastic_force[1, 2];
            Force[Tet[tet * 4 + 0]].z -= Elastic_force[2, 0] + Elastic_force[2, 1] + Elastic_force[2, 2];
        }

        Smooth_V(0.75f);

    	for(int i=0; i<number; i++)
    	{
    		//TODO: Update X and V here.
            V[i] = (V[i] + dt * Force[i] / mass) * damp;
            X[i] = X[i] + V[i] * dt;
            //TODO: (Particle) collision with floor.
            float d = Vector3.Dot(X[i] - P, N);
            if (d < 0.0f) // collision occur
            {
                float v_N_size = Vector3.Dot(V[i], N);
                // check velocity
                if (v_N_size < 0.0f)
                {
                    Vector3 v_N = v_N_size * N;
                    Vector3 v_T = V[i] - v_N;
                    Vector3 v_N_new = -1.0f * restitution * v_N;
                    float a = Math.Max(1.0f - friction * (1.0f + restitution) * v_N.magnitude / v_T.magnitude, 0.0f);
                    Vector3 v_T_new = a * v_T;
                    V[i] = v_N_new + v_T_new;
                }
            }

        }
    }

    // Update is called once per frame
    void Update()
    {
    	for(int l=0; l<10; l++)
    		 _Update();

    	// Dump the vertex array for rendering.
    	Vector3[] vertices = new Vector3[tet_number*12];
        int vertex_number=0;
        for(int tet=0; tet<tet_number; tet++)
        {
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        }
        Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.vertices  = vertices;
		mesh.RecalculateNormals ();
    }
}
