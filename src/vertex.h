#ifndef VERTEX_H
#define VERTEX_H

#include <QVector2D>
class Vertex
{
    float m_x,m_y;

    //int m_halfedgeIndex;

public:

    Vertex(float x = 0, float y = 0);
    Vertex(const Vertex& v);

    float x()const;
    float y()const;
    //int  halfedgeIndex()const;

    float& x();
    float& y();
    //int& halfedgeIndex();

    const Vertex& operator=(const Vertex&);
    QVector2D toVector2D();

};

#endif // VERTEX_H
