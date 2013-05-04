#ifndef CHEBUILDEROCTREEFROMPOINTCLOUD_H
#define CHEBUILDEROCTREEFROMPOINTCLOUD_H

#include "chebuilder.h"
#include <QVector>
#include <QPoint>
#include <QVector4D>

class OcTree
{
    class OcTreeNode
    {
    public:
        int xMin, xMax;
        int yMin, yMax;
        bool hasPoint;

        OcTreeNode *child00,*child01,*child11,*child10;

        OcTreeNode(int xm, int ym, int xM, int yM){
            xMin = xm; xMax = xM; yMin = ym; yMax = yM;
            child00 = child01 = child11 = child10 = (OcTreeNode*)NULL;
            hasPoint = false;
        }
        ~OcTreeNode(){
            if(child00)
                delete child00;
            if(child01)
                delete child01;
            if(child11)
                delete child11;
            if(child10)
                delete child10;
        }
        QVector<QPoint> getPoints(const QVector<QPoint>& points){
            QVector<QPoint> ret;
            for(int i = 0; i < points.size(); ++i)
            {
                if(points[i].x() >= xMin && points[i].x() < xMax &&
                        points[i].y() >= yMin && points[i].y() < yMax)
                    ret.append(points[i]);
            }
            return ret;
        }

    };


    OcTreeNode *root;

    float xMin, xMax;
    float yMin, yMax;

    int xMaxRes, yMaxRes;
    int subdivisionQuantityThreshold;

    QVector<Vertex> vertices;

public:
    typedef struct _TRI
    {
        int v[3];
    }TRI;

    OcTree(float _xMin, float _xMax, float _yMin, float _yMax);
    ~OcTree();

    void build(const QVector<QVector4D>& points);

    void getTriangles(QVector<Vertex>& vertices, QVector<TRI>& triangles);
    void recursivilyGetTriangles(OcTreeNode *node,QVector<Vertex>& vertices, QVector<TRI>& triangles, int **map);

private:
    void recursivilyBuild(OcTreeNode *node, QVector<QPoint> &points);


};

class CHEBuilderOctreeFromPointCloud : public CHEBuilder
{
    QVector<QVector4D> pointCloud;

    OcTree *tree;

public:
    CHEBuilderOctreeFromPointCloud(const QVector<QVector4D>& pointCloud);
    ~CHEBuilderOctreeFromPointCloud();

    virtual void build();
};

#endif // CHEBUILDEROCTREEFROMPOINTCLOUD_H
