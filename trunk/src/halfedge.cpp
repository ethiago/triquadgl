#include "halfedge.h"

HalfEdge::HalfEdge(int vertexIndex)
{
    m_vertexIndex = vertexIndex;
    m_twinIndex = -1;
}

HalfEdge::HalfEdge(const HalfEdge& h)
{
    m_vertexIndex = h.vertexIndex();
    m_twinIndex = h.twinIndex();
}

int  HalfEdge::vertexIndex()const
{
    return m_vertexIndex;
}

int& HalfEdge::vertexIndex()
{
    return m_vertexIndex;
}

int  HalfEdge::twinIndex()const
{
    return m_twinIndex;
}

int& HalfEdge::twinIndex()
{
    return m_twinIndex;
}

const HalfEdge& HalfEdge::operator=(const HalfEdge& h)
{
    m_vertexIndex = h.vertexIndex();
    m_twinIndex = h.twinIndex();

    return *this;
}

bool HalfEdge::hasTwin()const
{
    return m_twinIndex != -1;
}
