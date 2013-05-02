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

    int e0 = mesh.size();

    mesh.append(HalfEdge(idxV1)); //from idxV1 to idxV2
    mesh.append(HalfEdge(idxV2)); //from idxV2 to idxV3
    mesh.append(HalfEdge(idxV3)); //from idxV3 to idxV1

    configTwin(e0+0, idxV2);
    configTwin(e0+1, idxV3);
    configTwin(e0+2, idxV1);

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

void CompactHalfEdge::configTwin(int halfEdgeIdx, int destinyVertexIdx)
{
    for(int i = 0; i < mesh.size(); ++i)
    {
        if(mesh[i].vertexIndex() == destinyVertexIdx &&
                mesh[halfEdgeNext(i)].vertexIndex() == mesh[halfEdgeIdx].vertexIndex())
        {

            mesh[   i       ].twinIndex() = halfEdgeIdx;
            mesh[halfEdgeIdx].twinIndex() = i;

        }
    }
}
