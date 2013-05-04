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


    void addVertices(const QVector<Vertex>&);
    bool addTriangle(int idxV1, int idxV2, int idxV3);

    int halfEdgeNext(int heIdx);
    int halfEdgePrevious(int heIdx);

    int halfEdgeExternNext(int heIdx);
    int halfEdgeExternPrevious(int heIdx);

    int halfEdgeStarNext(int heIdx);
    int halfEdgeStarPrevious(int heIdx);
    int vertexId(int triangleId, int halfEdgeOffset);
    int sizeOfTriangles();
    int sizeOfVertices();
    int mostClosedVertex(const Vertex& v);

    const Vertex& vertex(int i)const;
    Vertex& vertex(int i);

    void addVertex(const Vertex& v);
    void joinVerticesAt(const Vertex& v);

    int nextExternHalfEdgeOf    (int vertexIdx);


private:
    void configTwin(int halfEdgeIdx, int destinyVertexIdx);

};

#endif // COMPACTHALFEDGE_H
