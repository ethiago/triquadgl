
#include "fitting.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>

#define MINIMUM 15

typedef qreal real;

real *points[2]; // coordinate x point
real b[3];
int n;

QVector<Quadric> processa(const QMatrix4x4& vInverse)
{
    gsl_matrix * A = gsl_matrix_alloc (n, 3*5);
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

    real chisq = 0.0;
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n,15);
    gsl_multifit_linear (A, B, x, cov, &chisq, work);
    gsl_multifit_linear_free (work);

    qDebug() << chisq;

    QVector<Quadric> resp;
    Quadric q;
    for(int i = 0; i < 3; ++i)
    {
        float x1 = gsl_vector_get(x, i*5);
        float x2 = gsl_vector_get(x, i*5 + 1);
        float x3 = gsl_vector_get(x, i*5 + 2);
        q.a_b_c = QVector3D(x1,x2,x3);
        x1 = gsl_vector_get(x, i*5 + 3);
        x2 = gsl_vector_get(x, i*5 + 4);
        q.d_e_f = QVector3D(x1, x2, -1.0);
        resp.push_back(q);
    }

    gsl_matrix_free (A);
    gsl_vector_free (B);
    gsl_vector_free (x);
    gsl_matrix_free (cov);

    return resp;
}

QVector<Quadric> processaQUAD()
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

    qDebug() << chisq;

    QVector<Quadric> resp;
    Quadric q;

    float x1 = gsl_vector_get(x, 0);
    float x2 = gsl_vector_get(x, 1);
    float x3 = gsl_vector_get(x, 2);
    q.a_b_c = QVector3D(x1,x2,x3);
    x1 = gsl_vector_get(x, 3);
    x2 = gsl_vector_get(x, 4);
    q.d_e_f = QVector3D(x1, x2, -1.0);
    resp.push_back(q);
    resp.push_back(q);
    resp.push_back(q);

    gsl_matrix_free (A);
    gsl_vector_free (B);
    gsl_vector_free (x);
    gsl_matrix_free (cov);

    return resp;
}

QVector<Quadric> fittingGSL(const QMatrix4x4& inv, const QVector<QVector2D>& inpoints)
{
    if(inpoints.size() < MINIMUM)
        return QVector<Quadric>();

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
