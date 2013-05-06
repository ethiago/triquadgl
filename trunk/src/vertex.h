#ifndef VERTEX_H
#define VERTEX_H

#include <QVector2D>
#include "quadric2d.h"

class Vertex
{
    float m_x,m_y;
    Quadric2D q;

    int m_halfedgeIndex;

public:

    Vertex(float x = 0, float y = 0);
    Vertex(const Vertex& v);
    Vertex(const QVector2D& v);

    float x()const;
    float y()const;
    int  halfedgeIndex()const;
    Quadric2D quadric()const;

    float& x();
    float& y();
    int& halfedgeIndex();
    Quadric2D& quadric();

    float distanceTo(const Vertex&) const;

    Vertex& operator=(const Vertex&);
    Vertex& operator=(const QVector2D&);

    QVector2D toVector2D()const;

    static bool projectVertexIntoSegment(const Vertex& _p, const Vertex& segmentV0, const Vertex& segmentV1, Vertex *ret);
    static float cross2D(const Vertex& p, const Vertex& p1, const Vertex&p2);
};



#endif // VERTEX_H
