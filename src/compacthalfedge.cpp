#include "compacthalfedge.h"

CompactHalfEdge::CompactHalfEdge()
{
}

CompactHalfEdge::CompactHalfEdge(const CompactHalfEdge& che)
{
    m_vertices = che.vertices();
    m_mesh = che.mesh();
}

CompactHalfEdge& CompactHalfEdge::operator=(const CompactHalfEdge& che)
{
    m_vertices = che.vertices();
    m_mesh = che.mesh();
    return *this;
}

void CompactHalfEdge::clear()
{
    m_vertices.clear();
    m_mesh.clear();
}

bool CompactHalfEdge::isEmpty()const
{
    return (m_vertices.size() == 0 || m_mesh.size() == 0);
}

const QVector<Vertex>&   CompactHalfEdge::vertices()const
{
    return m_vertices;
}
const QVector<HalfEdge>& CompactHalfEdge::mesh()const
{
    return m_mesh;
}

void CompactHalfEdge::addVertices(const QVector<Vertex>& vList)
{
    m_vertices = vList;
}

const Vertex& CompactHalfEdge::vertex(int i)const
{
    if(i < 0 || i >= m_vertices.size())
        return Vertex();

    return m_vertices[i];
}

Vertex& CompactHalfEdge::vertex(int i)
{
    Vertex v;
    if(i < 0 || i >= m_vertices.size())
        return v;

    return m_vertices[i];
}

bool CompactHalfEdge::addTriangle(int idxV1, int idxV2, int idxV3)
{
    if(idxV1 < 0 || idxV1 > m_vertices.size() ||
            idxV2 < 0 || idxV2 > m_vertices.size() ||
            idxV3 < 0 || idxV3 > m_vertices.size())
        return false;

    int e = m_mesh.size();

    m_mesh.append(HalfEdge(idxV1)); //from idxV1 to idxV2
    m_mesh.append(HalfEdge(idxV2)); //from idxV2 to idxV3
    m_mesh.append(HalfEdge(idxV3)); //from idxV3 to idxV1

    m_vertices[idxV1].halfedgeIndex() = e+0;
    m_vertices[idxV2].halfedgeIndex() = e+1;
    m_vertices[idxV3].halfedgeIndex() = e+2;

    configTwin(e+0, idxV2);
    configTwin(e+1, idxV3);
    configTwin(e+2, idxV1);

    return true;

}

int CompactHalfEdge::halfEdgeNext(int heIdx) const
{
    return 3*(heIdx/3) + (heIdx + 1)%3;
}

int CompactHalfEdge::halfEdgePrevious(int heIdx)const
{
    return 3*(heIdx/3) + (heIdx + 2)%3;
}

int CompactHalfEdge::halfEdgeExternNext(int heIdx) const
{
    if(m_mesh[heIdx].hasTwin())
        return -1;

    heIdx = halfEdgePrevious(heIdx);
    while(m_mesh[heIdx].hasTwin())
        heIdx = halfEdgePrevious(m_mesh[heIdx].twinIndex());

    return heIdx;
}

int CompactHalfEdge::halfEdgeExternPrevious(int heIdx)const
{
    if(m_mesh[heIdx].hasTwin())
        return -1;

    heIdx = halfEdgeNext(heIdx);
    while(m_mesh[heIdx].hasTwin())
        heIdx = halfEdgeNext(m_mesh[heIdx].twinIndex());

    return heIdx;

}

int CompactHalfEdge::halfEdgeStarNext(int heIdx) const
{
    return halfEdgeNext(m_mesh[heIdx].twinIndex());
}

int CompactHalfEdge::halfEdgeStarPrevious(int heIdx)const
{
    return m_mesh[halfEdgePrevious(heIdx)].twinIndex();
}

void CompactHalfEdge::configTwin(int halfEdgeIdx, int destinyVertexIdx)
{
    for(int i = 0; i < m_mesh.size(); ++i)
    {
        if(m_mesh[i].vertexIndex() == destinyVertexIdx &&
                m_mesh[halfEdgeNext(i)].vertexIndex() == m_mesh[halfEdgeIdx].vertexIndex())
        {

            m_mesh[   i       ].twinIndex() = halfEdgeIdx;
            m_mesh[halfEdgeIdx].twinIndex() = i;
            return;

        }
    }
}

int CompactHalfEdge::vertexId(int triangleId, int halfEdgeOffset) const
{
    return m_mesh[triangleId*3+halfEdgeOffset].vertexIndex();
}
int CompactHalfEdge::sizeOfTriangles()const
{
    return m_mesh.size()/3;
}

int CompactHalfEdge::sizeOfVertices() const
{
    return m_vertices.size();
}

void CompactHalfEdge::addVertex(const Vertex& p)
{
    int heIdx = -1;
    float dist = 0.0;
    for(int i = 0; i < m_mesh.size(); ++i)
    {
        Vertex v;
        if(!m_mesh[i].hasTwin() && Vertex::projectVertexIntoSegment(p, m_vertices[m_mesh[i].vertexIndex()],
                                                                  m_vertices[m_mesh[halfEdgeNext(i)].vertexIndex()], &v))
        {
            heIdx = i;
            dist = p.distanceTo(v);
            break;
        }
    }

    if(heIdx == -1)
        return;

    for(int i = heIdx+1; i < m_mesh.size(); ++i)
    {
        Vertex v;
        if(!m_mesh[i].hasTwin() && Vertex::projectVertexIntoSegment(p, m_vertices[m_mesh[i].vertexIndex()],
                                                                  m_vertices[m_mesh[halfEdgeNext(i)].vertexIndex()], &v))
        {
            if(dist > p.distanceTo(v))
            {
                heIdx = i;
                dist = p.distanceTo(v);
            }
        }
    }

    int vId = m_vertices.size();
    m_vertices.append(p);

    addTriangle(m_mesh[halfEdgeNext(heIdx)].vertexIndex(), m_mesh[heIdx].vertexIndex(), vId);

}

void CompactHalfEdge::joinVerticesAt(const Vertex& p)
{
    int vId = mostClosedVertex(p);
    if(vId < 0)
        return;
    int hId = nextExternHalfEdgeOf(vId);
    if (hId < 0)
        return;

    int hId1 = halfEdgeExternNext(halfEdgeExternNext(hId));
    int vId1 = m_mesh[hId1].vertexIndex();
    int hId2 = halfEdgeExternPrevious(halfEdgeExternPrevious(hId));
    int vId2 = m_mesh[hId2].vertexIndex();

    int vIdOther = vId2;
    float dist = p.distanceTo(m_vertices[vIdOther]);

    if(p.distanceTo(m_vertices[vId1]) < dist)
    {
        // change de values to mantain vId internaly before vIdOther
        vIdOther = vId;
        vId = vId1;
        hId = hId1;
    }

    int vIdMidle = m_mesh[halfEdgeExternPrevious(hId)].vertexIndex();

    addTriangle(vIdOther, vIdMidle, vId);

}

int CompactHalfEdge::mostClosedVertex(const Vertex& v) const
{
    if(m_vertices.isEmpty())
        return -1;

    int vId = 0;
    float mDist = v.distanceTo(m_vertices[vId]);
    for(int i = 1; i < m_vertices.size(); ++i)
    {
        if(v.distanceTo(m_vertices[i]) < mDist)
        {
            vId = i;
            mDist = v.distanceTo(m_vertices[i]);
        }
    }
    return vId;
}

int CompactHalfEdge::nextExternHalfEdgeOf(int vertexIdx)const
{
    int he = m_vertices[vertexIdx].halfedgeIndex();
    int aux = he;
    while(m_mesh[aux].hasTwin() && halfEdgeStarNext(aux) != he)
    {
        aux = halfEdgeStarNext(aux);
    }

    if(m_mesh[aux].hasTwin())
        return -1;
    return aux;
}

