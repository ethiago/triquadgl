#include "compacthalfedge.h"
#include <QVector4D>

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

int CompactHalfEdge::vertexNext(int hId)const
{
    return m_mesh[halfEdgeNext(hId)].vertexIndex();
}

int CompactHalfEdge::vertexPrevious(int hId) const
{
    return m_mesh[halfEdgePrevious(hId)].vertexIndex();
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
    if(!m_mesh[heIdx].hasTwin())
        return halfEdgeNext(halfEdgeExternNext(heIdx));

    return halfEdgeNext(m_mesh[heIdx].twinIndex());
}

int CompactHalfEdge::halfEdgeStarPrevious(int heIdx)const
{
    if(!m_mesh[halfEdgePrevious(heIdx)].hasTwin())
        return halfEdgeExternPrevious(halfEdgePrevious(heIdx));

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

void CompactHalfEdge::addVertexAsNewExternalTriangle(const Vertex& p)
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

QVector3D CompactHalfEdge::getBairicentricCoordinate(const Vertex& v, int triangleID)
{
    int v0 = vertexId(triangleID, 0);
    int v1 = vertexId(triangleID, 1);
    int v2 = vertexId(triangleID, 2);

    float t = Vertex::cross2D(m_vertices[v0], m_vertices[v1], m_vertices[v2]);
    float a0 = Vertex::cross2D(m_vertices[v1],m_vertices[v2],v);
    float a1 = Vertex::cross2D(m_vertices[v2],m_vertices[v0],v);
    float a2 = Vertex::cross2D(m_vertices[v0],m_vertices[v1],v);

    return QVector3D(a0/t,a1/t,a2/t);
}

void CompactHalfEdge::addVertexInExistingTriangle(const Vertex& v, int triangleID)
{

    if(!m_mesh[triangleID*3 + 0].hasTwin() || !m_mesh[triangleID*3 + 1].hasTwin() ||
            !m_mesh[triangleID*3 + 2].hasTwin())
        return;

    QVector3D bar = getBairicentricCoordinate(v, triangleID);
    int v0 = vertexId(triangleID,0);
    int v1 = vertexId(triangleID,1);
    int v2 = vertexId(triangleID,2);

    Vertex nv = m_vertices[v0]*bar.x() +
            m_vertices[v1]*bar.y() +
            m_vertices[v2]*bar.z();

    int inv = m_vertices.size();
    m_vertices.append(nv);
    addTriangle(v0, v1, inv);
    addTriangle(v1, v2, inv);
    addTriangle(v2, v0, inv);

    deleteTriangle(triangleID);

}

void CompactHalfEdge::addVertex(const Vertex& p)
{
    int tId = findTriangleWith(p);
    if(tId < 0)
        addVertexAsNewExternalTriangle(p);
    else
        addVertexInExistingTriangle(p,tId);
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

int CompactHalfEdge::findTriangleWith(const Vertex& v)
{
    for(int i = 0; i < sizeOfTriangles(); ++i)
    {
        Vertex v0, v1, v2;
        v0 = m_vertices[m_mesh[3*i+0].vertexIndex()];
        v1 = m_vertices[m_mesh[3*i+1].vertexIndex()];
        v2 = m_vertices[m_mesh[3*i+2].vertexIndex()];

        float a0 = Vertex::cross2D(v0,v1,v);
        float a1 = Vertex::cross2D(v1,v2,v);
        float a2 = Vertex::cross2D(v2,v0,v);

        if(a0*a1 > 0 && a0*a2 > 0)
            return i;

    }

    return -1;
}

void CompactHalfEdge::deleteTriangle(int tId)
{
    QVector<int> he;
    tId *=3;
    int vID = -1;

    for(int i = 0; i < 3 ; ++i)
    {
        if(m_mesh[tId+i].twinIndex() >= 0)
        {
            vID = m_mesh[halfEdgePrevious(tId+i)].vertexIndex();
            he.append( m_mesh[tId+i].twinIndex());
        }
    }

    if(he.size() == 0)
        return;


    if(he.size() > 1)
    {
        for(int i = 0; i < 3 ; ++i)
        {
            int vIdT = m_mesh[tId+i].vertexIndex();

            if(m_vertices[vIdT].halfedgeIndex() == tId || m_vertices[vIdT].halfedgeIndex() == tId+1 || m_vertices[vIdT].halfedgeIndex() == tId+2)
                m_vertices[vIdT].halfedgeIndex() = halfEdgeStarNext(m_vertices[vIdT].halfedgeIndex());
        }
    }else
    {
        int heb = m_vertices[vID].halfedgeIndex();
        int vIdT = vertexNext(heb);
        if(m_vertices[vIdT].halfedgeIndex() == tId || m_vertices[vIdT].halfedgeIndex() == tId+1 || m_vertices[vIdT].halfedgeIndex() == tId+2)
            m_vertices[vIdT].halfedgeIndex() = halfEdgeStarNext(m_vertices[vIdT].halfedgeIndex());

        vIdT = vertexPrevious(heb);
        if(m_vertices[vIdT].halfedgeIndex() == tId || m_vertices[vIdT].halfedgeIndex() == tId+1 || m_vertices[vIdT].halfedgeIndex() == tId+2)
            m_vertices[vIdT].halfedgeIndex() = halfEdgeStarNext(m_vertices[vIdT].halfedgeIndex());
    }


    for(int i = 0; i < he.size(); ++i)
    {
        m_mesh[he[i]].twinIndex() = -1;
    }
    m_mesh.erase(m_mesh.begin()+tId, m_mesh.begin()+tId+3);
    for(int i = 0; i < m_mesh.size(); ++i)
    {
        if(m_mesh[i].twinIndex() > tId)
            m_mesh[i].twinIndex() -= 3;
    }
    for(int i = 0; i < m_vertices.size(); ++i)
    {
        if(m_vertices[i].halfedgeIndex() > tId)
            m_vertices[i].halfedgeIndex() -= 3;
    }

    if(he.size() == 1)
    {
        m_vertices.erase(m_vertices.begin()+vID);
        for(int i = 0; i < m_mesh.size(); ++i)
        {
            if(m_mesh[i].vertexIndex() > vID)
                m_mesh[i].vertexIndex() -= 1;
        }
    }
}

QVector<HalfEdge> CompactHalfEdge::getAllInternalEdges()const
{
    QVector<HalfEdge> ret;
    QVector<bool> flags(m_mesh.size());

    for(int i = 0; i < flags.size(); ++i)
        flags[i] = false;


    for(int i = 0; i < m_mesh.size(); ++i)
    {
        if(m_mesh[i].hasTwin() && !flags[i])
        {
            ret.append(m_mesh[i]);
            flags[i] = true;
            flags[m_mesh[i].twinIndex()] = true;
        }
    }

    return ret;
}

int CompactHalfEdge::triangleIDFromHeId(int heid)const
{
    return heid%3;
}

void CompactHalfEdge::addToAllVertices(const QVector2D& v)
{
    for(int i = 0; i < m_vertices.size(); ++i)
    {
        m_vertices[i].x() += v.x();
        m_vertices[i].y() += v.y();
    }
}

QVector<QPair<Vertex, HalfEdge> > CompactHalfEdge::getIntersections(const QVector<QVector4D>& pontos)
{
    QVector<QPair<Vertex, HalfEdge> > ret;

    QVector<bool> flags(m_mesh.size());

    for(int i = 0; i < flags.size(); ++i)
        flags[i] = false;

    for(int i = 0; i < pontos.size() -1 ; ++i)
    {
        for(int j = 0; j < m_mesh.size(); ++j)
        {
            if(m_mesh[j].hasTwin() && !flags[j])
            {
                Vertex inter;
                Vertex Va = m_vertices[m_mesh[j].vertexIndex()];
                Vertex Vb = m_vertices[m_mesh[m_mesh[j].twinIndex() ].vertexIndex()];

                Vertex Pa = Vertex(pontos[i].toVector2D() );
                Vertex Pb = Vertex(pontos[i+1].toVector2D() );

                if(Vertex::intersection(Va,Vb,Pa,Pb, &inter))
                {
                    ret.append(qMakePair( inter, m_mesh[j] ));
                    flags[j] = true;
                    flags[m_mesh[j].twinIndex()] = true;
                }
            }
        }
        for(int i = 0; i < flags.size(); ++i)
            flags[i] = false;
    }

    return ret;
}

QVector2D CompactHalfEdge::getMaior()const
{
    return m_maior;
}

QVector2D CompactHalfEdge::getMenor()const
{
    return m_menor;
}

void CompactHalfEdge::maxMimCalc()
{
    if(m_vertices.size() <= 0)
        return;
    m_maior = m_vertices[0].toVector2D();
    m_menor = m_vertices[0].toVector2D();

    for(int i = 0; i < m_vertices.size(); ++i)
    {
        if(m_vertices[i].x() > m_maior.x())
            m_maior.setX(m_vertices[i].x());

        if(m_vertices[i].y() > m_maior.y())
            m_maior.setY(m_vertices[i].y());

        if(m_vertices[i].x() < m_menor.x())
            m_menor.setX(m_vertices[i].x());

        if(m_vertices[i].y() < m_menor.y())
            m_menor.setY(m_vertices[i].y());
    }
}

bool CompactHalfEdge::isVertexInternal(int vId)
{
    int heIdFirst = m_vertices[vId].halfedgeIndex();
    int heId = heIdFirst;
    do{
        if( heId == -1 || !m_mesh[heId].hasTwin())
            return false;

        heId = halfEdgeNext(m_mesh[heId].twinIndex());

    }while(heId != heIdFirst);

    return true;
}

QVector<int> CompactHalfEdge::getNeighbours(int vId)
{
    QVector<int> ret;

    int heIdFirst = m_vertices[vId].halfedgeIndex();
    int heId = heIdFirst;
    do{
        ret.append( m_mesh[m_mesh[heId].twinIndex()].vertexIndex());

        heId = halfEdgeNext(m_mesh[heId].twinIndex());

    }while(heId != heIdFirst);

    return ret;
}

QVector<int> CompactHalfEdge::getInternalVertices()
{
    QVector<int> ret;

    for(int i = 0; i < m_vertices.size(); ++i)
    {
        if(isVertexInternal(i))
            ret.append(i);
    }

    return ret;
}

int CompactHalfEdge::numberOfInternalVertices()
{
    return getInternalVertices().size();
}
