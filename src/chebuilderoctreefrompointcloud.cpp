#include "chebuilderoctreefrompointcloud.h"
#include <QDebug>


OcTree::OcTree(float _xMin, float _xMax, float _yMin, float _yMax)
{
    xMin = _xMin;
    xMax = _xMax;
    yMin = _yMin;
    yMax = _yMax;

    root = NULL;

    xMaxRes = 128;
    yMaxRes = 128;

    subdivisionQuantityThreshold = 40;
}

OcTree::~OcTree()
{
    if(root)
        delete root;
}

void OcTree::build(const QVector<QVector4D>& points)
{
    root = new OcTree::OcTreeNode(0, 0, xMaxRes, yMaxRes);
    QVector<QPoint> in;

    for(int i = 0; i < points.size(); ++i)
    {
        QVector2D v((points[i].x()-xMin)/(xMax-xMin), (points[i].y()-yMin)/(yMax-yMin));
        QPoint p(v.x()*(xMaxRes-1), v.y()*(yMaxRes-1));

        in.append(p);
    }

    recursivilyBuild(root, in);

    qDebug() << "####### FIM! ###########";

}

void OcTree::recursivilyBuild(OcTreeNode *node, QVector<QPoint>& points)
{
    qDebug() << points;

    qDebug() << "(" << node->xMin << "," << node->yMin << ") - " << "(" << node->xMax << "," << node->yMax << ")";

    if(points.size() < subdivisionQuantityThreshold)
    {
        if(points.size() > 0)
        {
            node->hasPoint = true;
            qDebug() << "Folha";
        }
        return;
    }

    int xMidle = (node->xMin+node->xMax)/2;
    int yMidle = (node->yMin+node->yMax)/2;

    node->child00 = new OcTreeNode(node->xMin, node->yMin, xMidle, yMidle );
    QVector<QPoint> p00 = node->child00->getPoints(points);

    node->child01 = new OcTreeNode(xMidle, node->yMin, node->xMax, yMidle );
    QVector<QPoint> p01 = node->child01->getPoints(points);

    node->child10 = new OcTreeNode(node->xMin, yMidle, xMidle, node->yMax );
    QVector<QPoint> p10 = node->child10->getPoints(points);

    node->child11 = new OcTreeNode(xMidle, yMidle, node->xMax, node->yMax );
    QVector<QPoint> p11 = node->child11->getPoints(points);

    if(p00.size()+p01.size()+p10.size()+p11.size() != points.size())
        qDebug() << "Deu Ruim!";

    points.clear();

    recursivilyBuild(node->child00, p00);
    recursivilyBuild(node->child01, p01);
    recursivilyBuild(node->child10, p10);
    recursivilyBuild(node->child11, p11);
}

void OcTree::getTriangles(QVector<Vertex>& vertices, QVector<TRI>& triangles)
{
    int **map = new int*[xMaxRes+1];
    for(int i = 0; i <= xMaxRes; ++i)
    {
        map[i] = new int[yMaxRes+1];
        for(int j = 0; j <= yMaxRes; ++j)
        {
            map[i][j] = -1;
        }
    }

    recursivilyGetTriangles(root, vertices, triangles, map);

    for(int i = 0; i <= xMaxRes; ++i)
    {
        delete[] map[i];
    }
    delete[] map;
}

void OcTree::recursivilyGetTriangles(OcTreeNode *node,QVector<Vertex>& vertices, QVector<TRI>& triangles, int **map)
{
    if(!node->child00)
    {
        if(!node->hasPoint)
            return;


        int vId00, vId01, vId10, vId11;
        int b = vertices.size();
        int i,j;

        i = node->xMin;  j = node->yMin;
        if(map[i][j] < 0)
        {
            Vertex v(xMin + i*( (xMax-xMin)/xMaxRes), yMin + j*( (yMax-yMin)/yMaxRes));
            vertices.append(v);
            map[i][j] = b++;
        }
        vId00 = map[i][j];

        i = node->xMax;  j = node->yMin;
        if(map[i][j] < 0)
        {
            Vertex v(xMin + i*( (xMax-xMin)/xMaxRes), yMin + j*( (yMax-yMin)/yMaxRes));
            vertices.append(v);
            map[i][j] = b++;
        }
        vId01 = map[i][j];

        i = node->xMax;  j = node->yMax;
        if(map[i][j] < 0)
        {
            Vertex v(xMin + i*( (xMax-xMin)/xMaxRes), yMin + j*( (yMax-yMin)/yMaxRes));
            vertices.append(v);
            map[i][j] = b++;
        }
        vId11 = map[i][j];

        i = node->xMin;  j = node->yMax;
        if(map[i][j] < 0)
        {
            Vertex v(xMin + i*( (xMax-xMin)/xMaxRes), yMin + j*( (yMax-yMin)/yMaxRes));
            vertices.append(v);
            map[i][j] = b++;
        }
        vId10 = map[i][j];

        TRI t;
        t.v[0] = vId00;
        t.v[1] = vId01;
        t.v[2] = vId11;
        triangles.append(t);
        t.v[0] = vId00;
        t.v[1] = vId11;
        t.v[2] = vId10;
        triangles.append(t);

        return;
    }

    recursivilyGetTriangles(node->child00, vertices, triangles, map);
    recursivilyGetTriangles(node->child01, vertices, triangles, map);
    recursivilyGetTriangles(node->child11, vertices, triangles, map);
    recursivilyGetTriangles(node->child10, vertices, triangles, map);
}

CHEBuilderOctreeFromPointCloud::CHEBuilderOctreeFromPointCloud(const QVector<QVector4D> & _pointCloud) : CHEBuilder(), pointCloud(_pointCloud)
{
    float xm, xM;
    float ym, yM;
    tree = (OcTree*)NULL;

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

    tree = new OcTree(xm, xM, ym, yM);
}

CHEBuilderOctreeFromPointCloud::~CHEBuilderOctreeFromPointCloud()
{
    if(tree)
        delete tree;
}

void CHEBuilderOctreeFromPointCloud::build()
{
    tree->build(pointCloud);

    QVector<Vertex> v;
    QVector<OcTree::TRI> t;

    tree->getTriangles(v,t);

    m_che.addVertices(v);
    for(int i = 0; i < t.size(); ++i)
        m_che.addTriangle(t[i].v[0], t[i].v[1], t[i].v[2]);
}
