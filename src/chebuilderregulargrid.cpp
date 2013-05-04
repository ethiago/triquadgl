#include "chebuilderregulargrid.h"

CHEBuilderRegularGrid::CHEBuilderRegularGrid(float _xMin, float _xMax, float _yMin, float _yMax, int _nSubX, int _nSubY) : CHEBuilder()
{
    xMin  = _xMin ;
    xMax  = _xMax ;
    yMin  = _yMin ;
    yMax  = _yMax ;
    nSubX = _nSubX;
    nSubY = _nSubY;
}

void CHEBuilderRegularGrid::build()
{
    m_che.clear();
    float dx = (xMax-xMin)/nSubX;
    float dy = (yMax-yMin)/nSubY;
    QVector<Vertex> list;

    for(int i = 0; i <= nSubX; ++i)
    {
        for(int j = 0; j <= nSubY; ++j)
        {
            list.append(Vertex(xMin + i*dx, yMin + j*dy));
        }
    }

    m_che.addVertices(list);

    for(int i = 0; i < nSubX; ++i)
    {
        for(int j = 0; j < nSubY; ++j)
        {
            int vIdx00 = (i+0)*(nSubY+1) + (j+0);
            int vIdx01 = (i+1)*(nSubY+1) + (j+0);
            int vIdx10 = (i+0)*(nSubY+1) + (j+1);
            int vIdx11 = (i+1)*(nSubY+1) + (j+1);

            m_che.addTriangle(vIdx00,vIdx01, vIdx11);
            m_che.addTriangle(vIdx00,vIdx11, vIdx10);
        }
    }
}
