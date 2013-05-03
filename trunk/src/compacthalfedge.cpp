#include "compacthalfedge.h"

CompactHalfEdge::CompactHalfEdge()
{
}

void CompactHalfEdge::addVertices(const QVector<Vertex>& vList)
{
    vertices = vList;
}

const Vertex& CompactHalfEdge::vertex(int i)const
{
    if(i < 0 || i >= vertices.size())
        return Vertex();

    return vertices[i];
}

Vertex& CompactHalfEdge::vertex(int i)
{
    Vertex v;
    if(i < 0 || i >= vertices.size())
        return v;

    return vertices[i];
}

bool CompactHalfEdge::addTriangle(int idxV1, int idxV2, int idxV3)
{
    if(idxV1 < 0 || idxV1 > vertices.size() ||
            idxV2 < 0 || idxV2 > vertices.size() ||
            idxV3 < 0 || idxV3 > vertices.size())
        return false;

    int e = mesh.size();

    mesh.append(HalfEdge(idxV1)); //from idxV1 to idxV2
    mesh.append(HalfEdge(idxV2)); //from idxV2 to idxV3
    mesh.append(HalfEdge(idxV3)); //from idxV3 to idxV1

    vertices[idxV1].halfedgeIndex() = e+0;
    vertices[idxV2].halfedgeIndex() = e+1;
    vertices[idxV3].halfedgeIndex() = e+2;

    configTwin(e+0, idxV2);
    configTwin(e+1, idxV3);
    configTwin(e+2, idxV1);

    return true;

}

int CompactHalfEdge::halfEdgeNext(int heIdx)
{
    return 3*(heIdx/3) + (heIdx + 1)%3;
}

int CompactHalfEdge::halfEdgePrevious(int heIdx)
{
    return 3*(heIdx/3) + (heIdx + 2)%3;
}

int CompactHalfEdge::halfEdgeExternNext(int heIdx)
{
    if(mesh[heIdx].hasTwin())
        return -1;

    heIdx = halfEdgePrevious(heIdx);
    while(mesh[heIdx].hasTwin())
        heIdx = halfEdgePrevious(mesh[heIdx].twinIndex());

    return heIdx;
}

int CompactHalfEdge::halfEdgeExternPrevious(int heIdx)
{
    if(mesh[heIdx].hasTwin())
        return -1;

    heIdx = halfEdgeNext(heIdx);
    while(mesh[heIdx].hasTwin())
        heIdx = halfEdgeNext(mesh[heIdx].twinIndex());

    return heIdx;

}

int CompactHalfEdge::halfEdgeStarNext(int heIdx)
{
    return halfEdgeNext(mesh[heIdx].twinIndex());
}

int CompactHalfEdge::halfEdgeStarPrevious(int heIdx)
{
    return mesh[halfEdgePrevious(heIdx)].twinIndex();
}

void CompactHalfEdge::configTwin(int halfEdgeIdx, int destinyVertexIdx)
{
    for(int i = 0; i < mesh.size(); ++i)
    {
        if(mesh[i].vertexIndex() == destinyVertexIdx &&
                mesh[halfEdgeNext(i)].vertexIndex() == mesh[halfEdgeIdx].vertexIndex())
        {

            mesh[   i       ].twinIndex() = halfEdgeIdx;
            mesh[halfEdgeIdx].twinIndex() = i;
            return;

        }
    }
}

int CompactHalfEdge::vertexId(int triangleId, int halfEdgeOffset)
{
    return mesh[triangleId*3+halfEdgeOffset].vertexIndex();
}
int CompactHalfEdge::sizeOfTriangles()
{
    return mesh.size()/3;
}

int CompactHalfEdge::sizeOfVertices()
{
    return vertices.size();
}

void CompactHalfEdge::addVertex(const Vertex& p)
{
    int heIdx = -1;
    float dist = 0.0;
    for(int i = 0; i < mesh.size(); ++i)
    {
        Vertex v;
        if(!mesh[i].hasTwin() && Vertex::projectVertexIntoSegment(p, vertices[mesh[i].vertexIndex()],
                                                                  vertices[mesh[halfEdgeNext(i)].vertexIndex()], &v))
        {
            heIdx = i;
            dist = p.distanceTo(v);
            break;
        }
    }

    if(heIdx == -1)
        return;

    for(int i = heIdx+1; i < mesh.size(); ++i)
    {
        Vertex v;
        if(!mesh[i].hasTwin() && Vertex::projectVertexIntoSegment(p, vertices[mesh[i].vertexIndex()],
                                                                  vertices[mesh[halfEdgeNext(i)].vertexIndex()], &v))
        {
            if(dist > p.distanceTo(v))
            {
                heIdx = i;
                dist = p.distanceTo(v);
            }
        }
    }

    int vId = vertices.size();
    vertices.append(p);

    addTriangle(mesh[halfEdgeNext(heIdx)].vertexIndex(), mesh[heIdx].vertexIndex(), vId);

}

void CompactHalfEdge::joinVerticesAt(const Vertex& p)
{
    int vId = mostClosedVertex(p);
    int hId = nextExternHalfEdgeOf(vId);
    if (hId < 0)
        return;

    int hId1 = halfEdgeExternNext(halfEdgeExternNext(hId));
    int vId1 = mesh[hId1].vertexIndex();
    int hId2 = halfEdgeExternPrevious(halfEdgeExternPrevious(hId));
    int vId2 = mesh[hId2].vertexIndex();

    int vIdOther = vId2;
    float dist = p.distanceTo(vertices[vIdOther]);

    if(p.distanceTo(vertices[vId1]) < dist)
    {
        // change de values to mantain vId internaly before vIdOther
        vIdOther = vId;
        vId = vId1;
        hId = hId1;
    }

    int vIdMidle = mesh[halfEdgeExternPrevious(hId)].vertexIndex();

    addTriangle(vIdOther, vIdMidle, vId);

}

int CompactHalfEdge::mostClosedVertex(const Vertex& v)
{
    int vId = 0;
    float mDist = v.distanceTo(vertices[vId]);
    for(int i = 1; i < vertices.size(); ++i)
    {
        if(v.distanceTo(vertices[i]) < mDist)
        {
            vId = i;
            mDist = v.distanceTo(vertices[i]);
        }
    }
    return vId;
}

int CompactHalfEdge::nextExternHalfEdgeOf(int vertexIdx)
{
    int he = vertices[vertexIdx].halfedgeIndex();
    int aux = he;
    while(mesh[aux].hasTwin() && halfEdgeStarNext(aux) != he)
    {
        aux = halfEdgeStarNext(aux);
    }

    if(mesh[aux].hasTwin())
        return -1;
    return aux;
}

//int CompactHalfEdge::previousExternHalfEdgeOf(int vertexIdx)
//{
//    int he = halfEdgePrevious(vertices[vertexIdx].halfedgeIndex());
//    int aux = he;
//    while(mesh[aux].hasTwin() && halfEdgePrevious(mesh[aux].twinIndex()))
//    {
//        halfEdgePrevious(mesh[aux].twinIndex());
//    }

//    if(mesh[aux].hasTwin())
//        return -1;
//    return aux;
//}
