#ifndef CHEBUILDERREGULARGRID_H
#define CHEBUILDERREGULARGRID_H

#include "chebuilder.h"

class CHEBuilderRegularGrid : public CHEBuilder
{
    float xMin  , xMax ;
    float yMin  , yMax ;
    int   nSubX , nSubY ;

public:
    CHEBuilderRegularGrid(float _xMin = -1.0, float _xMax = 1.0, float _yMin = -1.0, float _yMax = 1.0, int _nSubX = 10, int _nSubY = 10);

    virtual void build();
};



#endif // CHEBUILDERREGULARGRID_H
