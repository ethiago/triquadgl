#ifndef CHEBUILDERREGULARGRIDFROMPOINTCLOUD_H
#define CHEBUILDERREGULARGRIDFROMPOINTCLOUD_H

#include "chebuilder.h"
#include <QVector4D>
#include <QVector>
#include <QPoint>

class CHEBuilderRegularGridFromPointCloud : public CHEBuilder
{
    QVector<QVector4D> pointCloud;

    int nSubx;
    int nSuby;

public:
    CHEBuilderRegularGridFromPointCloud(const QVector<QVector4D>& _pointCloud, int _nSubx, int _nSubY);

    virtual void build();
};

#endif // CHEBUILDERREGULARGRIDFROMPOINTCLOUD_H
