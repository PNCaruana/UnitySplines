using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class nurbs_surface : MonoBehaviour
{
    public enum KNOT_TYPE
    {
        Normalized,
        Closed,
        Uniform
    }


    [SerializeField]
    public KNOT_TYPE KnotType;

    // Start is called before the first frame update

    /// <summary>
	/// number of control points along u axis (min 3)
	/// </summary>
    [SerializeField]
    public int uDim = 4;

    /// <summary>
	/// number of control points along v axis (min 3)
	/// </summary>
    [SerializeField]
    public int vDim = 4;

    /// <summary>
	/// Order in the u-direction
	/// </summary>
    [SerializeField]
    public int K = 3;

    /// <summary>
    /// Order in the v-direction
    /// </summary>
    [SerializeField]
    public int L = 3;

    private Vector3[,] controlPoints; // x,y as axis, with z as height
    private GameObject[,] controlObjects;

    

    void Start()
    {
        controlPoints = new Vector3[uDim , vDim];
        controlObjects = new GameObject[uDim, vDim];
        //init x,y with default z=0;
        //Once set, (x,y) are immutable and should not ever be modified
        int count = 0;
        for (int x=0; x<uDim; x++)
        {
            float u = (float)x / uDim;
            for (int y=0; y<vDim; y++)
            {
                float v = (float)y / vDim;

                Vector3 pt = new Vector3(u, 0, v);

                controlPoints[x, y] = pt;
                GameObject cp = GameObject.CreatePrimitive(PrimitiveType.Sphere);
                cp.transform.localScale = new Vector3(0.02f, 0.02f, 0.02f);
                cp.transform.position = pt;
                cp.transform.SetParent(transform);
                controlObjects[x, y] = cp;

                //Debug.Log("Created point at : " + (pt));

                count++;
            }
        }

        //Vector3 P1 = BS_Spline(0.05f);


    }

    //Test functions
    void testBS_Basis(float u, int i, int n, float[] U)
    {
        float result = BS_Basis(u, i, n, U);

        Debug.Log("Result (u = "+u+"): " + result);
    }


    string listToStr(int[] list)
    {
        string str = "<";
        foreach (int i in list)
        {
            str += i + ", ";
        }
        str.Remove(str.Length - 2);
        str += ">";

        return str;
    }

    string listToStr(float[] list)
    {
        string str = "<";
        foreach (float i in list)
        {
            str += i + ", ";
        }
        str.Remove(str.Length - 2);
        str += ">";

        return str;
    }

    // Update is called once per frame
    void Update()
    {
        float[] K = genKnotVector(7, 3, KNOT_TYPE.Closed);
        Debug.Log(listToStr(K));

        //Debug.Log(BS_Basis(0.9999f, 6, 3, K));
        drawSurface();
    }

    void drawSurface() {
        for (int i = 0; i < uDim; i++)
        {
            for (int j = 0; j < vDim; j++)
            {
                controlPoints[i, j] = controlObjects[i, j].transform.position; //update control points to the point entities
            }

        }

        float u = 0;
        float v = 0;
        float t = 0;

        int N = 20; //how many segments

        float dt = (float)1 / N;


        for (int i = 0; i < N; i++)
        {
            //Vector3 P1 = BS_Spline(t);
            //Vector3 P2 = BS_Spline(t+dt);
            //Debug.DrawLine(P1, P2, Color.red, 0.01f);

            for (int j = 0; j < N; j++)
            {
                u = i * dt;// + 0.0001f;
                v = j * dt;// + 0.0001f;

                //Debug.Log("(u,v) " + u + ", " + v);

                if (i == 0) u = 0;
                if (i == N) u = 1;
                if (j == 0) v = 0;
                if (j == N) v = 1;

                float u2 = u + dt;
                float v2 = v + dt;

                if (i == N - 1)
                {
                    u2 = 0.9999f;
                }

                if (j == N - 1)
                {
                    v2 = 0.9999f;
                }

                //yield return new WaitForSeconds(0.001f);

                Vector3 P1 = BS_Surface(u, v);
                Vector3 P2 = BS_Surface(u2, v);
                Vector3 P3 = BS_Surface(u, v2);

                Debug.DrawLine(P1, P2, Color.red, 0.01f);
                Debug.DrawLine(P1, P3, Color.red, 0.01f);
            }

        }


        //draw end edges
        for (int i = 0; i < N; i++)
        {
            u = i * dt;
            v = 0.9999f;

            float u2 = u + dt;

            if (i == N - 1)
            {
                u2 = 0.9999f;
            }

            Vector3 P1 = BS_Surface(u, v);
            Vector3 P2 = BS_Surface(u2, v);
            Debug.DrawLine(P1, P2, Color.red, 0.01f);
        }

        for (int j = 0; j < N; j++)
        {
            u = 0.9999f;
            v = j * dt;

            float v2 = v + dt;

            if (j == N - 1)
            {
                v2 = 0.9999f;
            }

            Vector3 P1 = BS_Surface(u, v);
            Vector3 P2 = BS_Surface(u, v2);
            Debug.DrawLine(P1, P2, Color.red, 0.01f);
        }
    }

    
    //Interface Functions
    public void setHeights(float[,] heights) {
        if (heights.Length != controlPoints.Length) {
            Debug.LogWarning("!! Incorrect size of float[,] heights specified ("+heights.Length+"), expected length should be # of control points");

        }

        for (int i=0; i<uDim; i++)
        {
            for (int j = 0; j < vDim; j++) {
                Vector3 oldPt = controlPoints[i, j];
                controlObjects[i, j].transform.position = new Vector3(oldPt.x, -heights[uDim -1 - i, j], oldPt.z);
            }
        }
    }
    //==============================================================================

    Vector3 bezierSurface(float u, float v) {
        Vector3 P = Vector3.zero;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                P += controlPoints[i, j] * blend(u, i, 3) * blend(v, j, 3);

            }
        }

        return P;
    }

    void drawBezier(Vector3 start, Vector3 p0, Vector3 p1, Vector3 end) {
        float t = 0;
        float dt = 0.01f;
        for (int i = 0; i < 100; i++) {
            t = i * dt;
            Vector3 A = cubicBezier(start, p0, p1, end, t);
            Vector3 B = cubicBezier(start, p0, p1, end, t +dt);
            Debug.DrawLine(A, B, Color.red, 0.05f);
        }
    }

    //parametric value <t>, control point index <i> and degree <n>
    //Bezier Blend Function
    float blend(float t, int i, int n) {
        float bin = nCr(n, i);

        return bin * Mathf.Pow(t, i) * Mathf.Pow(1 - t, n - i);
    }

    Vector3 cubicBezier(Vector3 start,  Vector3 p0, Vector3 p1, Vector3 end, float t) {
        Vector3 A = quadBezier(start, p1, p0, t);
        Vector3 B = quadBezier(p0, end, p1, t);
        Vector3 AB = Vector3.Lerp(A, B, t);

        return AB;
    }

    Vector3 quadBezier(Vector3 x0, Vector3 x1, Vector3 p0, float t) {
        Vector3 A = Vector3.Lerp(x0, p0, t);
        Vector3 B = Vector3.Lerp(p0, x1, t);
        Vector3 AB = Vector3.Lerp(A, B, t);

        return AB;
    }

    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // BSplines ===================================================================================
    // VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

    /*
     * parameter <u>
     * basis number <i> 
     * Knot Vector <U>
     * Degree <n>
     */
    float BS_Basis(float u, int i,  int n, float[] U) {

        float basis;
        // End of recursive basis functions for degree 0
        // N<i,0>(u)
        if (n == 0)
        {
            

            // degree 0 basis is simply a step function if <t> within knot
            if ((U[i] <= u) && (u < U[i + 1]))
            {
                basis = 1;
            }
            else
            {
                basis = 0;
            }
        }
        else {
            // N<i,n>(u)

            //calculate w<i,n,u>
            

            float lc = weight(i, n, u, U); //Left coefficient
            float rc = 1 - weight(i+1, n, u, U); //Right coefficient

            float bL = BS_Basis(u, i, n - 1, U);
            float bR = BS_Basis(u, i+1, n - 1, U);

            float lhs = lc * bL;
            float rhs = rc * bR;

            //Sometimes lc and rc are NAN, but we dont care cause when that happens basis should be 0 anyways

            /*
            if (bL == 0)
            {
                lhs = 0;
            }
            else {
                lhs = lc * bL;
            }

            if (bR == 0)
            {
                rhs = 0;
            }
            else
            {
                rhs = rc * bR;
            }
            */

            basis = lhs + rhs;
        }

        if (u == 0.05f) {
            //Debug.Log("Basis("+i+", "+n+") = " + basis);
        }
        

        return basis;
    }

    //coeficient for basis recursion
    float weight(int i, int n, float u, float[] U) {

        float w = 0;

        if (U[i + n] != U[i])
        {
            w = (u - U[i]) / (U[i + n] - U[i]);
        }

        return w;
    }

    //First derivative of Basis function
    float BS_Basis_Prime(float u, int i, int n, float[] U)
    {
        return alpha(i, n, U) * BS_Basis(u, i, n - 1, U) - alpha(i + 1, n, U) * BS_Basis(u, i + 1, n - 1, U);
    }

    //Alpha coefficient for derivative
    //point i, degree n
    float alpha(int i, int n, float[] U) {
        if (U[i + n] == U[i]) {
            return 0;
        }
        else
        {
            return n / (U[i + n] - U[i]);
        }
    }

    public Vector3 surfaceGradient(float u, float v) {
        // dz/dx = dz/du / dx/du + dz/dv / dx/dv
        float[] U = genKnotVector(uDim, K, KnotType);
        float[] V = genKnotVector(vDim, L, KnotType);
        //Compute the jacobian
        Vector3 du = Vector3.zero; //du = <dx/du, dy/du, dz/du>
        Vector3 dv = Vector3.zero; //db = <dx/dv, dy/dv, dz/dv>

        for (int i = 0; i < uDim; i++) {
            for (int j =0; j < vDim; j++)
            {
                du += controlPoints[i, j] * BS_Basis_Prime(u, i, K, U) * BS_Basis(v, j, L, V);
            }
        }

        
        for (int i = 0; i < uDim; i++)
        {
            for (int j = 0; j < vDim; j++)
            {
                dv += controlPoints[i, j] * BS_Basis(u, i, K, U) * BS_Basis_Prime(v, j, L, V);
            }
        }

        Debug.Log("du " + du);
        Debug.Log("dv " + dv);
        //Now we want to compute dy/dx and dy/dz

        float dydx;
        float dydz;

        dydx = du.y / du.x; //chain rule
        dydz = dv.y / dv.z;

        return new Vector3(dydx, 0, dydz);
    }

    float[] genKnotVector(int numPts, int degree, KNOT_TYPE K) {

        int h = numPts + degree;
        float[] U = new float[h+1];
        for (int i = 0; i < h+1; i++) {

            if (K == KNOT_TYPE.Normalized)
            {
                U[i] = i * (1f / h);
            }
            else if (K == KNOT_TYPE.Uniform)
            {
                U[i] = i;
            }
            else if (K == KNOT_TYPE.Closed)
            {
                if (i+1 <= degree) {
                    U[i] = 0;
                } else if (i+1 >= (numPts+1)) {
                    U[i] = 1;
                }
                else {
                    U[i] = (float)(i - degree) / (numPts - degree);
                }

                
            }
        }

        return U;

    }

    Vector3 BS_Surface(float u, float v) {
        Vector3 P = Vector3.zero;
        float[] U = genKnotVector(uDim, K, KnotType);
        float[] V = genKnotVector(vDim, L, KnotType);

        for (int i = 0; i < uDim; i++)
        {
            for (int j = 0; j < vDim; j++)
            {
                P += controlPoints[i, j] * BS_Basis(u, i, K, U) * BS_Basis(v, j, L, V);

            }
        }
        if (P == Vector3.zero)
        {
            //Debug.LogWarning("Zero Vec @ u,v = " + u + ", " + v);
        }

        return P;
    }

    Vector3 BS_Spline(float t)
    {

        Vector3 P = Vector3.zero;
        float[] U = genKnotVector(uDim, K, KnotType);
        Debug.Log(listToStr(U));
        for (int i = 0; i < uDim; i++)
        {
            P += controlPoints[i, 0] * BS_Basis(t, i, K, U);
            //Debug.Log("i = " + i + ", P = " + P);
        }

        if (P == Vector3.zero) {
            //Debug.LogWarning("Zero Vec @ t = " + t);
        }

        return P;

    }

    float nCr(int n, int r) {
        int facN = factorial(n);
        int facR = factorial(r);
        int facNR = factorial(n - r);

        return facN / (facR * facNR);
    }

    int factorial(int n) {
        if (n == 0)
        {
            return 1;
        }

        return n * factorial(n - 1);
    }
}
