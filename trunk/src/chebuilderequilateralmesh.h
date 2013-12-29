#ifndef CHEBUILDEREQUILATERALMESH_H
#define CHEBUILDEREQUILATERALMESH_H

#include "chebuilder.h"
#include <QVector>
#include <QVector4D>

class CHEBuilderEquilateralMesh : public CHEBuilder
{

    float ymin, xmin, xmax;
    int nSubX, nSubY;
    float l,h;

    int index(int i, int j);

public:
    CHEBuilderEquilateralMesh(const QVector<QVector4D>& _pointCloud, int _nSubx);

     virtual void build();
};

#endif // CHEBUILDEREQUILATERALMESH_H
