#ifndef COMPACTHALFEDGE_H
#define COMPACTHALFEDGE_H

#include "vertex.h"
#include "halfedge.h"
#include <QVector>

class CompactHalfEdge
{
    QVector<Vertex> m_vertices;
    QVector<HalfEdge> m_mesh;

public:
    CompactHalfEdge();
    CompactHalfEdge(const CompactHalfEdge& che);

    CompactHalfEdge& operator=(const CompactHalfEdge& che);


    const QVector<Vertex>&   vertices()const;
    const QVector<HalfEdge>& mesh()const;

    void clear();
    bool isEmpty()const;


    void addVertices(const QVector<Vertex>&);
    bool addTriangle(int idxV1, int idxV2, int idxV3);

    int halfEdgeNext(int heIdx)const;
    int halfEdgePrevious(int heIdx) const;

    int halfEdgeExternNext(int heIdx)const;
    int halfEdgeExternPrevious(int heIdx) const;

    int halfEdgeStarNext(int heIdx)const;
    int halfEdgeStarPrevious(int heIdx) const;
    int vertexId(int triangleId, int halfEdgeOffset)const;
    int sizeOfTriangles() const;
    int sizeOfVertices()const;
    int mostClosedVertex(const Vertex& v)const;

    const Vertex& vertex(int i)const;
    Vertex& vertex(int i);

    void addVertex(const Vertex& v);
    void joinVerticesAt(const Vertex& v);

    int nextExternHalfEdgeOf    (int vertexIdx) const;


private:
    void configTwin(int halfEdgeIdx, int destinyVertexIdx);

};

#endif // COMPACTHALFEDGE_H
