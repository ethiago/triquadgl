#include "compacthalfedge.h"

CompactHalfEdge::CompactHalfEdge()
{
}

void CompactHalfEdge::addVertices(const QVector<Vertex>& vList)
{
    vertices = vList;
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

void CompactHalfEdge::addVertex(const Vertex& v)
{
    float mDist = v.distanceTo(vertices[0]);
    int mDistIdx = 0;
    for(int i = 1; i < vertices.size(); ++i)
    {
        if( v.distanceTo(vertices[i]) < mDist)
        {
            mDist = v.distanceTo(vertices[i]);
            mDistIdx = i;
        }
    }

    int nehe = nextExternHalfEdgeOf(mDistIdx);
    int pehe = previousExternHalfEdgeOf(mDistIdx);


    int newVertexId = vertices.size();
    vertices.append(v);
    if(v.distanceTo(mesh[halfEdgeNext(nehe)].vertexIndex()) < v.distanceTo(mesh[pehe].vertexIndex()) )
        addTriangle(mesh[halfEdgeNext(nehe)].vertexIndex(), mesh[nehe].vertexIndex(), newVertexId);
    else
        addTriangle(mesh[halfEdgeNext(pehe)].vertexIndex(), mesh[pehe].vertexIndex(), newVertexId);

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

int CompactHalfEdge::previousExternHalfEdgeOf(int vertexIdx)
{
    int he = halfEdgePrevious(vertices[vertexIdx].halfedgeIndex());
    int aux = he;
    while(mesh[aux].hasTwin() && halfEdgePrevious(mesh[aux].twinIndex()))
    {
        halfEdgePrevious(mesh[aux].twinIndex());
    }

    if(mesh[aux].hasTwin())
        return -1;
    return aux;
}
