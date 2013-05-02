#ifndef HALFEDGE_H
#define HALFEDGE_H

class HalfEdge
{
    int m_vertexIndex;

    int m_twinIndex;

public:
    HalfEdge(int vertexIndex = -1);
    HalfEdge(const HalfEdge&);

    int  vertexIndex()const;
    int& vertexIndex();

    int  twinIndex()const;
    int& twinIndex();

    bool hasTwin()const;

    const HalfEdge& operator=(const HalfEdge&);

};

#endif // HALFEDGE_H
