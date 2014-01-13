#include "chebuilderequilateralmesh.h"
#include "qmath.h"

CHEBuilderEquilateralMesh::CHEBuilderEquilateralMesh(const QVector<QVector4D> &pointCloud, int _nSubx) :
    nSubX(_nSubx), nSubY(0)
{
    if(pointCloud.empty())
        return;


    float ymax = ymin = pointCloud[0].y();
    xmin = xmax = pointCloud[0].x();


    for(int i = 1; i < pointCloud.size(); ++i)
    {
        if(pointCloud[i].y() < ymin)
            ymin = pointCloud[i].y();

        if(pointCloud[i].y() > ymax)
            ymax = pointCloud[i].y();

        if(pointCloud[i].x() < xmin)
            xmin = pointCloud[i].x();

        if(pointCloud[i].x() > xmax)
            xmax = pointCloud[i].x();
    }

    float hX = xmax-xmin;

    xmax += hX*0.005;
    xmin -= hX*0.005;


    l = (xmax - xmin)/nSubX;
    h = (sqrt(3.0)*l)/2.0;

    nSubY = qCeil((ymax - ymin)/h);

    ymin -= (nSubY*h - (ymax - ymin))/2.0;
}

int CHEBuilderEquilateralMesh::index(int x, int y)
{
    return y*(nSubX+1) + x + (y/2);
}

void CHEBuilderEquilateralMesh::build()
{
    if(nSubY == 0)
        return;

    QVector<Vertex> vertices;

    float y = ymin;
    for(int j = 0; j <= nSubY ; ++j)
    {
        float x = j%2==0 ? xmin : xmin - l/2.0;
        for(int i = 0; i <= nSubX + (j%2); ++i)
        {
            vertices.append(Vertex(x,y));
            x+=l;
        }
        y += h;
    }

    m_che.addVertices(vertices);

    for(int j = 0; j < nSubY ; ++j)
    {
        for(int i = 0; i < nSubX ; ++i)
        {
            int idx0 = index(i,j);
            int idxNX = index(i+1,j);
            int idxNY = index(i,j+1);
            int idxNXY = index(i+1,j+1);
            if(j%2 == 0)
            {
                m_che.addTriangle(idx0, idxNXY, idxNY);
                m_che.addTriangle(idx0, idxNX, idxNXY);
            }
            else
            {
                m_che.addTriangle(idx0,idxNX, idxNY);
                m_che.addTriangle(idxNX,idxNXY, idxNY);
            }
        }
        int i = nSubX;
        int idx0 = index(i,j);
        int idxNX = index(i+1,j);
        int idxNY = index(i,j+1);
        int idxNXY = index(i+1,j+1);
        if(j%2 == 0)
        {
            m_che.addTriangle(idx0, idxNXY, idxNY);
        }
        else
        {
            m_che.addTriangle(idx0,idxNX, idxNY);
        }
    }

}
