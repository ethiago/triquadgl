#include "quadric2d.h"

Quadric2D::Quadric2D()
{
    a = 0.0;
    b = 0.0;
    c = 0.0;
    d = 0.0;
    e = 0.0;
    f = 0.0;
}


Quadric2D::Quadric2D(float x2, float xy, float x, float y2, float y, float k)
{
    a = x2;
    b = xy/2.0;
    c = x/2.0;
    d = y2;
    e = y/2.0;
    f = k;
}

QMatrix3x3 Quadric2D::toMatrixForm()
{
    QMatrix3x3 m;
    m(0,0) = a;
    m(0,1) = m(1,0) = b;
    m(0,2) = m(2,0) = c;
    m(1,1) = d;
    m(1,2) = m(2,1) = e;
    m(2,2) = f;
    return m;
}
