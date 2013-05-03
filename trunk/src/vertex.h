#ifndef VERTEX_H
#define VERTEX_H

#include <QVector2D>
class Vertex
{
    float m_x,m_y;

    int m_halfedgeIndex;

public:

    Vertex(float x = 0, float y = 0);
    Vertex(const Vertex& v);

    float x()const;
    float y()const;
    int  halfedgeIndex()const;

    float& x();
    float& y();
    int& halfedgeIndex();

    float distanceTo(const Vertex&) const;

    const Vertex& operator=(const Vertex&);
    QVector2D toVector2D();

    static bool projectVertexIntoSegment(const Vertex& _p, const Vertex& segmentV0, const Vertex& segmentV1, Vertex *ret);

};

#endif // VERTEX_H
