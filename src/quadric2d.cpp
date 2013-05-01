#include "quadric2d.h"

Quadric2D::Quadric2D()
{
    m_a = 0.0;
    m_b = 0.0;
    m_c = 0.0;
    m_d = 0.0;
    m_e = 0.0;
    m_f = 0.0;
}


Quadric2D::Quadric2D(float x2, float xy, float x, float y2, float y, float k)
{
    m_a = x2;
    m_b = xy/2.0;
    m_c = x/2.0;
    m_d = y2;
    m_e = y/2.0;
    m_f = k;
}

QMatrix3x3 Quadric2D::toMatrixForm() const
{
    QMatrix3x3 m;
    m(0,0) = m_a;
    m(0,1) = m(1,0) = m_b;
    m(0,2) = m(2,0) = m_c;
    m(1,1) = m_d;
    m(1,2) = m(2,1) = m_e;
    m(2,2) = m_f;
    return m;
}

float Quadric2D::a() const
{
    return m_a;
}

float Quadric2D::b() const
{
    return m_b;
}

float Quadric2D::c() const
{
    return m_c;
}

float Quadric2D::d() const
{
    return m_d;
}

float Quadric2D::e() const
{
    return m_e;
}

float Quadric2D::f() const
{
    return m_f;
}

QVector3D Quadric2D::abc() const
{
    return QVector3D(m_a,m_b,m_c);
}

QVector3D Quadric2D::def()const
{
    return QVector3D(m_d,m_e,m_f);
}
