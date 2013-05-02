#ifndef COMPACTHALFEDGE_H
#define COMPACTHALFEDGE_H

#include "vertex.h"
#include "halfedge.h"
#include <QVector>

class CompactHalfEdge
{
    QVector<Vertex> vertices;
    QVector<HalfEdge> mesh;

public:
    CompactHalfEdge();

    void addVertices(const QVector<Vertex>&);
    bool addTriangle(int idxV1, int idxV2, int idxV3);

    int halfEdgeNext(int heIdx);
    int halfEdgePrevious(int heIdx);

private:
    void configTwin(int halfEdgeIdx, int destinyVertexIdx);

};

#endif // COMPACTHALFEDGE_H
