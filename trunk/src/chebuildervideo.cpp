#include "chebuildervideo.h"

CHEBuilderVideo::CHEBuilderVideo(int type) : CHEBuilder()
{
    m_type = type;
}

void CHEBuilderVideo::build()
{
    m_che.clear();
    QVector<Vertex> hev;
    Vertex v;
    switch(m_type)
    {
    case 0:
        v = Vertex( -2, -1);
        hev.append(v);
        v = Vertex(  2, -1);
        hev.append(v);
        v = Vertex( 0, 1);
        hev.append(v);
        break;
    case 1:
        v = Vertex( -2, -1);
        v.quadric() = Quadric2D(1,0,0,1,0,-0.3);
        hev.append(v);
        v = Vertex(  2, -1);
        v.quadric() = Quadric2D(1,0,0,1,0,-0.3);
        hev.append(v);
        v = Vertex( 0, 1);
        v.quadric() = Quadric2D(1,0,0,1,0,-0.3);
        hev.append(v);
        break;
    case 2:
        v = Vertex( -2, -1);
        v.quadric() = Quadric2D(2,0,0,0,1,0);
        hev.append(v);
        v = Vertex(  2, -1);
        v.quadric() = Quadric2D(2,0,0,0,1,0);
        hev.append(v);
        v = Vertex( 0, 1);
        v.quadric() = Quadric2D(2,0,0,0,1,0);
        hev.append(v);
        break;
    case 3:
        v = Vertex( -2, -1);
        v.quadric() = Quadric2D(2,0,0,0,1,0);
        hev.append(v);
        v = Vertex(  2, -1);
        v.quadric() = Quadric2D(2,0,0,0,1,0);
        hev.append(v);
        v = Vertex( 0, 1);
        v.quadric() = Quadric2D(1,0,0,1,0,-0.3);
        hev.append(v);
        break;
    }

    m_che.addVertices(hev);

    m_che.addTriangle(0,1,2);
}
