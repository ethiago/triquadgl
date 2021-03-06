#include "vertex.h"
#include "qmath.h"
#include <QMatrix4x4>

Vertex::Vertex(float x , float y )
{
    m_x = x;
    m_y = y;
    m_halfedgeIndex = -1;
    q = ZERO;
}

Vertex::Vertex(const QVector2D& v)
{
    m_x = v.x();
    m_y = v.y();
    m_halfedgeIndex = -1;
    q = ZERO;
}

Vertex::Vertex(const Vertex& v)
{
    m_x = v.x();
    m_y = v.y();
    m_halfedgeIndex = v.halfedgeIndex();
    q = v.quadric();
}

float Vertex::x()const
{
    return m_x;
}
float& Vertex::x()
{
    return m_x;
}

float Vertex::y()const
{
    return m_y;
}
float& Vertex::y()
{
    return m_y;
}

int  Vertex::halfedgeIndex()const
{
    return m_halfedgeIndex;
}
int&  Vertex::halfedgeIndex()
{
    return m_halfedgeIndex;
}

Quadric2D Vertex::quadric()const
{
    return q;
}

Quadric2D& Vertex::quadric()
{
    return q;
}

Vertex& Vertex::operator=(const Vertex& v)
{
    m_x = v.x();
    m_y = v.y();
    m_halfedgeIndex = v.halfedgeIndex();
    q = v.quadric();

    return *this;
}

Vertex& Vertex::operator=(const QVector2D& v)
{
    m_x = v.x();
    m_y = v.y();

    return *this;
}

Vertex Vertex::operator+(const Vertex& v)const
{
    Vertex nv(m_x+v.x(), m_y+v.y());
    nv.quadric() = quadric() + v.quadric();
    return nv;
}

Vertex Vertex::operator*(float f) const
{
    Vertex v(m_x * f, m_y * f);
    v.quadric() = quadric()*f;
    return v;
}

QVector2D Vertex::toVector2D()const
{
    return QVector2D(m_x,m_y);
}

float Vertex::distanceTo(const Vertex& v)const
{
    return sqrt( (v.x() - m_x)*(v.x() - m_x) + (v.y() - m_y)*(v.y() - m_y));
}

bool Vertex::projectVertexIntoSegment(const Vertex& _p, const Vertex& segmentV0, const Vertex& segmentV1, Vertex* ret)
{
    QVector2D p(_p.x(), _p.y());
    QVector2D v0(segmentV0.x(), segmentV0.y());
    QVector2D v1(segmentV1.x(), segmentV1.y());
    QVector2D d(v1-v0);

    float t = QVector2D::dotProduct(p-v0, d)/QVector2D::dotProduct(d, d);

    QVector2D r = v0 + t*d;
    ret->x() = r.x();
    ret->y() = r.y();

    if(t >= 0.0 && t <= 1.0)
        return true;
    else
        return false;
}

float Vertex::cross2D(const Vertex& p, const Vertex& p1, const Vertex& p2)
{
    float v1X = p1.x()-p.x();
    float v1Y = p1.y()-p.y();
    float v2X = p2.x()-p.x();
    float v2Y = p2.y()-p.y();

    return v1X*v2Y - v2X*v1Y;
}

bool Vertex::intersection(const Vertex& Va, const Vertex& Vb, const Vertex& Pa, const Vertex& Pb, Vertex *inter)
{
    QMatrix4x4 M;

    M(0,0) = Vb.x() - Va.x(); M(0,1) = Pa.x() - Pb.x();
    M(1,0) = Vb.y() - Va.y(); M(1,1) = Pa.y() - Pb.y();

    bool i;
    M = M.inverted(&i);

    if(!i)
        return false;

    QVector4D B(Pa.x() - Va.x(), Pa.y() - Va.y(), 0, 0);
    B = M*B;

    if( 0 <= B.x() && B.x() <= 1.0 && 0 <= B.y() && B.y() <= 1.0)
    {
        *inter = Va.toVector2D() + B.x()*( Vb.toVector2D() - Va.toVector2D() );
        return true;
    }

    return false;
}
