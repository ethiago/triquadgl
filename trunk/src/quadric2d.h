#ifndef QUADRIC2D_H
#define QUADRIC2D_H

#include <QMatrix3x3>
#include <QVector3D>

class Quadric2D
{
    float m_a,m_b,m_c,m_d,m_e,m_f;

public:
    Quadric2D();
    Quadric2D(float x2, float xy, float x, float y2, float y, float k);

    QMatrix3x3 toMatrixForm()const;

    float a()const;
    float b()const;
    float c()const;
    float d()const;
    float e()const;
    float f() const;

    QVector3D abc()const;
    QVector3D def() const;


};

#endif // QUADRIC2D_H
