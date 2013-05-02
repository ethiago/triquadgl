#ifndef VERTEX_H
#define VERTEX_H


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

    const Vertex& operator=(const Vertex&);

};

#endif // VERTEX_H
