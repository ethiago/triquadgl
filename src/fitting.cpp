
#include "fitting.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>

#define MINIMUM 15

typedef qreal real;

real *points[2]; // coordinate x point
real b[3];
int n;

QVector<Quadric2D> processa(const QMatrix4x4& vInverse)
{
    gsl_matrix * A = gsl_matrix_calloc (n, 3*5);
    real bary[3];

    for (int i = 0; i < n; ++i)
    {
        QVector4D b(points[0][i], points[1][i], 1.0, 1.0);
        b = vInverse * b;
        bary[0] = b.x();
        bary[1] = b.y();
        bary[2] = b.z();
        for (int j = 0; j < 3; ++j)
        {
            gsl_matrix_set (A, i, j*5    ,     points[0][i]*points[0][i]*bary[j]); //x^2
            gsl_matrix_set (A, i, j*5 + 1, 2.0*points[0][i]*points[1][i]*bary[j]); //2xy
            gsl_matrix_set (A, i, j*5 + 2, 2.0*points[0][i]             *bary[j]); //2x
            gsl_matrix_set (A, i, j*5 + 3,     points[1][i]*points[1][i]*bary[j]); //y^2
            gsl_matrix_set (A, i, j*5 + 4, 2.0*points[1][i]             *bary[j]); //2y
        }
    }

    gsl_matrix *cov = gsl_matrix_alloc (15, 15);
    gsl_vector * B = gsl_vector_alloc(n);
    gsl_vector * x = gsl_vector_alloc(15);
    gsl_vector_set_all(B, 1.0);

    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n,15);
    real chisq = 0.0;
//    gsl_multifit_linear (A, B, x, cov, &chisq, work);
    real tol = 0.0001 ;
    size_t rank = 0 ;
    gsl_multifit_linear_svd( A,B, tol, &rank, x, cov, &chisq, work);
    gsl_multifit_linear_free (work);

    //qDebug() << chisq;
    //qDebug() << rank ;

    QVector<Quadric2D> resp;
    for(int i = 0; i < 3; ++i)
    {
        float x1 = gsl_vector_get(x, i*5);
        float x2 = gsl_vector_get(x, i*5 + 1);
        float x3 = gsl_vector_get(x, i*5 + 2);
        float x4 = gsl_vector_get(x, i*5 + 3);
        float x5 = gsl_vector_get(x, i*5 + 4);

        resp.push_back(Quadric2D(x1,x2,x3,x4,x5,-1.0));
    }

    gsl_matrix_free (A);
    gsl_vector_free (B);
    gsl_vector_free (x);
    gsl_matrix_free (cov);

    return resp;
}

QVector<Quadric2D> processaQUAD()
{
    gsl_matrix * A = gsl_matrix_alloc (n, 5);

    for (int i = 0; i < n; ++i)
    {
        gsl_matrix_set (A, i, 0,     points[0][i]*points[0][i]); //x^2
        gsl_matrix_set (A, i, 1, 2.0*points[0][i]*points[1][i]); //2xy
        gsl_matrix_set (A, i, 2, 2.0*points[0][i]             ); //2x
        gsl_matrix_set (A, i, 3,     points[1][i]*points[1][i]); //y^2
        gsl_matrix_set (A, i, 4, 2.0*points[1][i]             ); //2y
    }

    gsl_matrix *cov = gsl_matrix_alloc (5, 5);
    gsl_vector * B = gsl_vector_alloc(n);
    gsl_vector * x = gsl_vector_alloc(5);
    gsl_vector_set_all(B, 1.0);

    real chisq = 0.0;
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n,5);
    gsl_multifit_linear (A, B, x, cov, &chisq, work);
    gsl_multifit_linear_free (work);

    //qDebug() << chisq;

    QVector<Quadric2D> resp;

    float x1 = gsl_vector_get(x, 0);
    float x2 = gsl_vector_get(x, 1);
    float x3 = gsl_vector_get(x, 2);
    float x4 = gsl_vector_get(x, 3);
    float x5 = gsl_vector_get(x, 4);

    resp.push_back(Quadric2D(x1,x2,x3,x4,x5,-1.0));
    resp.push_back(Quadric2D(x1,x2,x3,x4,x5,-1.0));
    resp.push_back(Quadric2D(x1,x2,x3,x4,x5,-1.0));

    gsl_matrix_free (A);
    gsl_vector_free (B);
    gsl_vector_free (x);
    gsl_matrix_free (cov);

    return resp;
}

QVector<Quadric2D> fittingGLOBAL(gsl_matrix * A, gsl_vector * B, float f)
{

    gsl_matrix *cov = gsl_matrix_alloc (A->size2, A->size2);
    gsl_vector * x = gsl_vector_alloc(A->size2);

    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (A->size1,A->size2);
    real chisq = 0.0;
    real tol = 0.000001 ;
    size_t rank = 0 ;
    gsl_multifit_linear_svd( A,B, tol, &rank, x, cov, &chisq, work);
    //gsl_multifit_linear(A,B, x, cov, &chisq, work);
    gsl_multifit_linear_free (work);

    //qDebug() << chisq;
    //qDebug() << rank ;

    QVector<Quadric2D> resp;

    for(unsigned int i = 0; i < A->size2/5; ++i)
    {
        float x1 = gsl_vector_get(x, i*5);
        float x2 = gsl_vector_get(x, i*5 + 1); // b/2
        float x3 = gsl_vector_get(x, i*5 + 2); // c/2
        float x4 = gsl_vector_get(x, i*5 + 3);
        float x5 = gsl_vector_get(x, i*5 + 4); // e/2

        resp.push_back(Quadric2D(x1,x2,x3,x4,x5,f));
    }

    gsl_vector_free (x);
    gsl_matrix_free (cov);

    return resp;
}

 gsl_vector* simpleFitting(gsl_matrix * A, gsl_vector * B)
{
    gsl_matrix *cov = gsl_matrix_alloc (A->size2, A->size2);
    gsl_vector * x = gsl_vector_alloc(A->size2);

    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (A->size1,A->size2);
    real chisq = 0.0;
    real tol = 0.000001 ;
    size_t rank = 0 ;
    gsl_multifit_linear_svd( A,B, tol, &rank, x, cov, &chisq, work);
    qDebug() << "rank: " << rank << ", chisq: " << chisq;
    //gsl_multifit_linear(A,B, x, cov, &chisq, work);
    gsl_multifit_linear_free (work);

    return x;
}

QVector<Quadric2D> fittingGLOBAL_flivre(gsl_matrix * A, gsl_vector * B)
{

    gsl_matrix *cov = gsl_matrix_alloc (A->size2, A->size2);
    gsl_vector * x = gsl_vector_alloc(A->size2);

    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (A->size1,A->size2);
    real chisq = 0.0;
    real tol = 0.000001 ;
    size_t rank = 0 ;
    gsl_multifit_linear_svd( A,B, tol, &rank, x, cov, &chisq, work);
    //gsl_multifit_linear(A,B, x, cov, &chisq, work);
    gsl_multifit_linear_free (work);

    //qDebug() << chisq;
    //qDebug() << rank ;

    QVector<Quadric2D> resp;
    for(unsigned int i = 0; i < A->size2/6; ++i)
    {
        float x1 = gsl_vector_get(x, i*6);
        float x2 = gsl_vector_get(x, i*6 + 1);// b/2
        float x3 = gsl_vector_get(x, i*6 + 2);// c/2
        float x4 = gsl_vector_get(x, i*6 + 3);
        float x5 = gsl_vector_get(x, i*6 + 4);// e/2
        float x6 = gsl_vector_get(x, i*6 + 5);

        resp.push_back(Quadric2D(x1,x2,x3,x4,x5,x6));
    }

    gsl_vector_free (x);
    gsl_matrix_free (cov);

    return resp;
}

QVector<Quadric2D> fittingGSL(const QMatrix4x4& inv, const QVector<QVector2D>& inpoints)
{
    if(inpoints.size() < MINIMUM)
    {
        QVector<Quadric2D> qs ;
        qs.push_back( ZERO ) ;
        qs.push_back( ZERO ) ;
        qs.push_back( ZERO ) ;
        return qs;
    }

    n = inpoints.size();
    points[0] = new real[n];
    points[1] = new real[n];

    for(int i = 0; i < n; ++i)
    {
        points[0][i] = inpoints[i].x();
        points[1][i] = inpoints[i].y();
    }

    return processa(inv);
}
