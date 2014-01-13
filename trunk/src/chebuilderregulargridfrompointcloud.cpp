#include "chebuilderregulargridfrompointcloud.h"

CHEBuilderRegularGridFromPointCloud::CHEBuilderRegularGridFromPointCloud(const QVector<QVector4D>& _pointCloud, int _nSubx, int _nSubY) :
    CHEBuilder(), pointCloud(_pointCloud), nSubx(_nSubx), nSuby(_nSubY)
{
}

void CHEBuilderRegularGridFromPointCloud::build()
{
    float xm, xM;
    float ym, yM;

    if(pointCloud.isEmpty())
        return;

    xm = xM = pointCloud[0].x();
    ym = yM = pointCloud[0].y();

    for(int i = 1; i < pointCloud.size(); ++i)
    {
        if(pointCloud[i].x() < xm)
            xm = pointCloud[i].x();
        if(pointCloud[i].y() < ym)
            ym = pointCloud[i].y();

        if(pointCloud[i].x() > xM)
            xM = pointCloud[i].x();
        if(pointCloud[i].y() > yM)
            yM = pointCloud[i].y();
    }

    float lx = (xM - xm);
    float ly = (yM - ym);

    float perc = 0.05;

    xm -= lx*0.005;
    xM += lx*0.005;
    ym -= ly*perc;
    yM += ly*perc;

    bool **buckets = new bool*[nSubx];
    for(int i = 0; i < nSubx; ++i)
    {
        buckets[i] = new bool[nSuby];
        for(int j = 0; j < nSuby; ++j)
            buckets[i][j] = false;
    }

    int **map = new int*[nSubx+1];
    for(int i = 0; i <= nSubx; ++i)
    {
        map[i] = new int[nSuby+1];
        for(int j = 0; j <= nSuby; ++j)
            map[i][j] = -1;
    }

    for(int i = 0; i < pointCloud.size(); ++i)
    {
        int xIdx = nSubx*(pointCloud[i].x() - xm)/(xM-xm);
        int yIdx = nSuby*(pointCloud[i].y() - ym)/(yM-ym);
        buckets[xIdx][yIdx] = true;
    }

    QVector<Vertex> vertices;
    QVector<int> tri;
    int vId = 0;
    for(int i = 0; i < nSubx; ++i)
    {
        for(int j = 0; j < nSuby; ++j)
        {
            if(buckets[i][j])
            {
                int vId00, vId01, vId10, vId11;
                int pi,pj;

                pi = i;   pj = j;
                if(map[pi][pj] < 0)
                {
                    Vertex v(xm + pi*( (xM-xm)/nSubx), ym + pj*( (yM-ym)/nSuby));
                    vertices.append(v);
                    map[pi][pj] = vId++;
                }
                vId00 = map[pi][pj];

                pi = i+1;   pj = j;
                if(map[pi][pj] < 0)
                {
                    Vertex v(xm + pi*( (xM-xm)/nSubx), ym + pj*( (yM-ym)/nSuby));
                    vertices.append(v);
                    map[pi][pj] = vId++;
                }
                vId01 = map[pi][pj];

                pi = i+1;   pj = j+1;
                if(map[pi][pj] < 0)
                {
                    Vertex v(xm + pi*( (xM-xm)/nSubx), ym + pj*( (yM-ym)/nSuby));
                    vertices.append(v);
                    map[pi][pj] = vId++;
                }
                vId11 = map[pi][pj];

                pi = i;   pj = j+1;
                if(map[pi][pj] < 0)
                {
                    Vertex v(xm + pi*( (xM-xm)/nSubx), ym + pj*( (yM-ym)/nSuby));
                    vertices.append(v);
                    map[pi][pj] = vId++;
                }
                vId10 = map[pi][pj];

                tri.append(vId00); tri.append(vId01); tri.append(vId11);
                tri.append(vId00); tri.append(vId11); tri.append(vId10);
            }
        }
    }


    for(int i = 0; i < nSubx; ++i)
        delete[] buckets[i];
    delete[] buckets;

    for(int i = 0; i <= nSubx; ++i)
        delete[] map[i];
    delete[] map;

    m_che.addVertices(vertices);
    for(int i = 0; i < tri.size()/3; ++i)
        m_che.addTriangle(tri[i*3+0], tri[i*3+1], tri[i*3+2]);

}
