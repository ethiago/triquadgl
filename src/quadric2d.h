#ifndef QUADRIC2D_H
#define QUADRIC2D_H

#include <QMatrix3x3>

class Quadric2D
{
    float a,b,c,d,e,f;

public:
    Quadric2D();
    Quadric2D(float x2, float xy, float x, float y2, float y, float k);

    QMatrix3x3 toMatrixForm();
};

#endif // QUADRIC2D_H
