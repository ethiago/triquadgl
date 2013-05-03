#include "vertex.h"
#include "qmath.h"

Vertex::Vertex(float x , float y )
{
    m_x = x;
    m_y = y;
}

Vertex::Vertex(const Vertex& v)
{
    m_x = v.x();
    m_y = v.y();
    m_halfedgeIndex = v.halfedgeIndex();
}

float Vertex::x()const
{
    return m_x;
}

float Vertex::y()const
{
    return m_y;
}

int  Vertex::halfedgeIndex()const
{
    return m_halfedgeIndex;
}

float& Vertex::x()
{
    return m_x;
}

float& Vertex::y()
{
    return m_y;
}

int&  Vertex::halfedgeIndex()
{
    return m_halfedgeIndex;
}

const Vertex& Vertex::operator=(const Vertex& v)
{
    m_x = v.x();
    m_y = v.y();
    m_halfedgeIndex = v.halfedgeIndex();

    return *this;
}
QVector2D Vertex::toVector2D()
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
