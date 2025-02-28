/**
 * general math routines
 * @author Tobias Weber <tweber@ill.fr>
 * @date 29-jun-19
 * @license see 'LICENSE' file
 */

public class Calc
{
    // ------------------------------------------------------------------------
    /**
     * r_i = x_i + y_i
     */
    /*@
        requires x.length == y.length;
        ensures \result.length == x.length;
        ensures (\forall int i; 
                 0<=i && i<\result.length; 
                 \result[i] - (x[i] + y[i]) < 1e-6 ||
                 -(\result[i] - (x[i] + y[i])) < 1e-6);
    @*/
    public static double[] add(double[] x, double[] y)
        throws Exception
    {
        final int dim = x.length;
        if(dim != y.length)
            throw new Exception("Vector sizes are not compatible.");

        double[] r = new double[dim];
        for(int i=0; i<dim; ++i)
            r[i] = x[i] + y[i];

        return r;
    }


    /**
     * r_i = x_i - y_i
     */
    /*@
        requires x.length == y.length;
        ensures \result.length == x.length;
        ensures (\forall int i;
                 0<=i && i<\result.length;
                 \result[i] - (x[i] - y[i]) < 1e-6 ||
                 -(\result[i] - (x[i] - y[i])) < 1e-6);
    @*/
    public static double[] sub(double[] x, double[] y)
        throws Exception
    {
        final int dim = x.length;
        if(dim != y.length)
            throw new Exception("Vector sizes are not compatible.");

        double[] r = new double[dim];
        for(int i=0; i<dim; ++i)
            r[i] = x[i] - y[i];

        return r;
    }


    /**
     * r_i = x_i * d
     */
    /*@
        ensures \result.length == x.length;
        ensures (\forall int i; 
                 0<=i && i<\result.length;
                 \result[i] - x[i]*d < 1e-6 ||
                 -(\result[i] - x[i]*d) < 1e-6);
    @*/
    public static double[] mul(double[] x, double d)
        throws Exception
    {
        final int dim = x.length;

        double[] r = new double[dim];
        for(int i=0; i<dim; ++i)
            r[i] = x[i] * d;

        return r;
    }


    /**
     * r_i = d * x_i
     */
    /*@
        ensures \result.length == x.length;
        ensures (\forall int i; 
                 0<=i && i<\result.length;
                 \result[i] - x[i]*d < 1e-6 ||
                 -(\result[i] - x[i]*d) < 1e-6);
    @*/
    public static double[] mul(double d, double[] x)
        throws Exception
    {
        return mul(x, d);
    }


    /**
     * r_i = x_i / d
     */
    /*@
        ensures \result.length == x.length;
        ensures (\forall int i; 
                 0<=i && i<\result.length;
                 \result[i] - x[i]/d < 1e-6 ||
                 -(\result[i] - x[i]/d) < 1e-6);
    @*/
    public static double[] div(double[] x, double d)
        throws Exception
    {
        return mul(x, 1./d);
    }


    /**
     * R_ij = A_ij * d
     */
    public static double[][] mul(double[][] A, double d)
        throws Exception
    {
        final int dim1 = A.length;
        final int dim2 = A[0].length;

        double[][] R = new double[dim1][dim2];
        for(int i=0; i<dim1; ++i)
            for(int j=0; j<dim2; ++j)
                R[i][j] = A[i][j] * d;

        return R;
    }


    /**
     * R_ij = d * A_ij
     */
    public static double[][] mul(double d, double[][] A)
        throws Exception
    {
        return mul(A, d);
    }


    /**
     * R_ij = A_ij / d
     */
    public static double[][] div(double[][] A, double d)
        throws Exception
    {
        return mul(A, 1./d);
    }


    /**
     * d = x_i y_i
     */
    public static double dot(double[] x, double[] y)
        throws Exception
    {
        final int dim = x.length;
        if(dim != y.length)
            throw new Exception("Vector sizes are not compatible.");

        double d = 0.;
        for(int i=0; i<dim; ++i)
            d += x[i]*y[i];

        return d;
    }


    /**
     * r_i = eps_ijk x_j y_k
     */
    public static double[] cross(double[] x, double[] y)
        throws Exception
    {
        final int dim = x.length;
        if(dim != y.length)
            throw new Exception("Vector sizes are not compatible.");
        if(dim != 3)
            throw new Exception("Dimension has to be 3.");

        double[] d = new double[]
        {
            x[1]*y[2] - x[2]*y[1],
            x[2]*y[0] - x[0]*y[2],
            x[0]*y[1] - x[1]*y[0]
        };

        return d;
    }


   /**
     * length of vector x
     */
    public static double norm_2(double[] x)
        throws Exception
    {
        return Math.sqrt(dot(x, x));
    }


    /**
     * d_i = M_ij x_j
     */
    public static double[] dot(double[][] M, double[] x)
        throws Exception
    {
        final int dim1 = M.length;
        final int dim2 = M[0].length;
        if(dim2 != x.length)
            throw new Exception("Matrix and vector sizes are not compatible.");

        double[] d = new double[dim1];
        for(int i=0; i<dim1; ++i)
            for(int j=0; j<dim2; ++j)
                d[i] += M[i][j]*x[j];

        return d;
    }


    /**
     * R_ik = A_ij x B_jk
     */
    public static double[][] dot(double[][] A, double[][] B)
        throws Exception
    {
        final int dim1 = A.length;
        final int dim2 = B[0].length;
        final int diminner = A[0].length;
        if(diminner != B.length)
            throw new Exception("Matrix sizes are not compatible.");

        double[][] R = new double[dim1][dim2];
        for(int i=0; i<dim1; ++i)
            for(int k=0; k<dim2; ++k)
                for(int j=0; j<diminner; ++j)
                    R[i][k] += A[i][j] * B[j][k];

        return R;
    }


    /**
     * R_ij = A_ij + B_ij
     */
    public static double[][] add(double[][] A, double[][] B)
        throws Exception
    {
        final int dim1 = A.length;
        final int dim2 = A[0].length;
        if(dim1 != B.length || dim2 != B[0].length)
            throw new Exception("Matrix sizes are not compatible.");

        double[][] R = new double[dim1][dim2];
        for(int i=0; i<dim1; ++i)
            for(int j=0; j<dim2; ++j)
                R[i][j] = A[i][j] + B[i][j];

        return R;
    }


    /**
     * R_ij = A_ij - B_ij
     */
    public static double[][] sub(double[][] A, double[][] B)
        throws Exception
    {
        final int dim1 = A.length;
        final int dim2 = A[0].length;
        if(dim1 != B.length || dim2 != B[0].length)
            throw new Exception("Matrix sizes are not compatible.");

        double[][] R = new double[dim1][dim2];
        for(int i=0; i<dim1; ++i)
            for(int j=0; j<dim2; ++j)
                R[i][j] = A[i][j] - B[i][j];

        return R;
    }


    /**
     * M_ij -> M_ji
     */
    public static double[][] transpose(double[][] M)
        throws Exception
    {
        final int dim1 = M.length;
        final int dim2 = M[0].length;

        double[][] R = new double[dim2][dim1];
        for(int i=0; i<dim1; ++i)
            for(int j=0; j<dim2; ++j)
                R[j][i] = M[i][j];

        return R;
    }


    /**
     * submatrix
     * @param M matrix
     * @param i row to delete
     * @param j column to delete
     * @return submatrix
     */
    public static double[][] submat(double[][] M, int i, int j)
    {
        final int dim1 = M.length;
        final int dim2 = M[0].length;

        double[][] R = new double[dim1-1][dim2-1];

        int _i2 = 0;
        for(int _i=0; _i<dim1; ++_i)
        {
            if(_i == i)
                continue;

            int _j2 = 0;
            for(int _j=0; _j<dim2; ++_j)
            {
                if(_j == j)
                    continue;

                R[_i2][_j2] = M[_i][_j];
                ++_j2;
            }
            ++_i2;
        }

        return R;
    }


    public static double[][] zero(int dim)
    {
        double[][] M = new double[dim][dim];

        for(int i=0; i<dim; ++i)
            for(int j=0; j<dim; ++j)
                M[i][j] = 0.;

        return M;
    }


    public static double[][] ident(int dim)
    {
        double[][] M = zero(dim);

        for(int i=0; i<dim; ++i)
            M[i][i] = 1.;

        return M;
    }


    /**
     * determinant
     * @param M matrix
     * @return determinant
     * @throws Exception
     */
    public static double det(double[][] M)
        throws Exception
    {
        final int dim1 = M.length;
        final int dim2 = M[0].length;

        if(dim1 != dim2)
            throw new Exception("Expecting a square matrix.");

        if(dim1 <= 0)
            return 0.;
        else if(dim1 == 1)
            return M[0][0];
        else if(dim1 == 2)
            return M[0][0]*M[1][1] - M[0][1]*M[1][0];

        double d = 0.;
        int i = 0;
        for(int j=0; j<dim2; ++j)
        {
            if(Math.abs(M[i][j]) < Double.MIN_VALUE)
                continue;

            double sgn = (((i+j) % 2) == 0) ? 1. : -1.;
            d += sgn * M[i][j] * det(submat(M, i,j));
        }

        return d;
    }


    /**
     * inverse
     * @param M matrix
     * @return inverse
     * @throws Exception
     */
    public static double[][] inv(double[][] M)
        throws Exception
    {
        final int dim1 = M.length;
        final int dim2 = M[0].length;

        if(dim1 != dim2)
            throw new Exception("Expecting a square matrix.");

        double d = det(M);
        double[][] I = new double[dim2][dim1];

        for(int i=0; i<dim1; ++i)
        {
            for(int j=0; j<dim2; ++j)
            {
                double sgn = ((i+j) % 2) == 0 ? 1. : -1.;
                I[j][i] = sgn * det(submat(M, i,j)) / d;
            }
        }

        return I;
    }


    /**
     * rotates a vector around an axis using Rodrigues' formula
     * see: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
     */
    public static double[] rotate(double[] _axis, double[] vec, double phi)
        throws Exception
    {
        double[] axis = div(_axis, norm_2(_axis));

        double s = Math.sin(phi);
        double c = Math.cos(phi);

        double[] vec1 = mul(vec, c);
        double[] vec2 = mul(mul(axis, dot(vec, axis)), (1.-c));
        double[] vec3 = mul(cross(axis, vec), s);

        return add(add(vec1, vec2), vec3);
    }


    /**
     * dot product in fractional coordinates
     */
    public static double dot(double[] x, double[] y, double[][] metric)
        throws Exception
    {
        return dot(x, dot(metric, y));
    }


    /**
     * levi-civita in cartesian coordinates
     */
    public static double levi(int i, int j, int k)
        throws Exception
    {
        int cols[] = new int[] {i, j, k};
        double[][] B = ident(3);
        double[][] M = ident(3);

        for(int idx1=0; idx1<B.length; ++idx1)
            for(int idx2=0; idx2<B[0].length; ++idx2)
                M[idx1][idx2] = B[idx1][cols[idx2]];

        return det(M);
    }


    /**
     * levi-civita symbol in fractional coordinates
     */
    public static double levi(int i, int j, int k, double[][] B)
        throws Exception
    {
        final int dim1 = B.length;
        final int dim2 = B[0].length;

        int cols[] = new int[] {i, j, k};
        double[][] M = new double[dim1][dim2];

        for(int idx1=0; idx1<dim1; ++idx1)
            for(int idx2=0; idx2<dim2; ++idx2)
                M[idx1][idx2] = B[idx1][cols[idx2]];

        return det(M);
    }


    /**
     * cross product in fractional coordinates
     */
    public static double[] cross(double[] a, double[] b, double[][] B)
        throws Exception
    {
        final int dim = B.length;

        double[][] metric_inv = inv(get_metric(B));
        double[] vec = new double[dim];

        for(int l=0; l<dim; ++l)
            for(int i=0; i<dim; ++i)
                for(int j=0; j<dim; ++j)
                    for(int k=0; k<dim; ++k)
                        vec[l] = levi(i,j,k, B) * a[j] * b[k] * metric_inv[l][i];
        return vec;
    }


    /**
     * angle between peaks in fractional coordinates
     */
    public static double angle(double[] x, double[] y, double[][] metric)
        throws Exception
    {
        double len_x = Math.sqrt(dot(x, x, metric));
        double len_y = Math.sqrt(dot(y, y, metric));

        double c = dot(x, y, metric) / (len_x * len_y);
        return Math.acos(c);
    }


    /**
     * get metric from crystal B matrix
     * basis vectors are in the columns of B, i.e. the second index
     */
    public static double[][] get_metric(double[][] B)
        throws Exception
    {
        return dot(transpose(B), B);
    }
    // ------------------------------------------------------------------------


    /*
     * testing openjml:
     * static checking: java -jar /opt/openjml/openjml.jar -esc Calc.java
     * dynamic checking: java -jar /opt/openjml/openjml.jar -rac Calc.java
     * java -cp .:/opt/openjml/jmlruntime.jar Calc
     */
    /*public static void main(String[] args)
    {
        try
        {
            double[] x = new double[]{1., 2., 3.};
            double[] y = new double[]{3., 4., 5.};
            double[] z1 = Calc.add(x, y);
            double[] z2 = Calc.sub(x, y);
        }
        catch(Exception ex)
        {
            System.err.println(ex);
        }
    }*/
}
