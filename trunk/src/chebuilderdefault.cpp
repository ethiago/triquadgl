#include "chebuilderdefault.h"

CHEBuilderDefault::CHEBuilderDefault() : CHEBuilder()
{
}

void CHEBuilderDefault::build()
{
    m_che.clear();

    QVector<Vertex> hev;
    hev.append(Vertex( -1.0, -1.0));
    hev.append(Vertex(  1.0, -1.0));
    hev.append(Vertex(  1.0,  1.0));
    hev.append(Vertex( -1.0,  1.0));
    hev.append(Vertex(  2.0,  0.0));
    hev.append(Vertex( -2.0,  0.0));

    m_che.addVertices(hev);

    m_che.addTriangle(0,1,2);
    m_che.addTriangle(2,3,0);
    m_che.addTriangle(4,2,1);
    m_che.addTriangle(5,0,3);
}
