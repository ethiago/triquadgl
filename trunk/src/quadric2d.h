#ifndef QUADRIC2D_H
#define QUADRIC2D_H

#include <QMatrix3x3>
#include <QVector3D>

#define ZERO     Quadric2D(0.0, 0.0, 0.0, 0.0, 0.0, 1.0)

class Quadric2D
{
    float m_a,m_b,m_c,m_d,m_e,m_f;

public:
    Quadric2D();
    Quadric2D(float , float , float, float , float,  float);
    Quadric2D(const Quadric2D&);


    QMatrix3x3 toMatrixForm()const;

    float a()const;
    float b()const;
    float c()const;
    float d()const;
    float e()const;
    float f()const;

    float& a();
    float& b();
    float& c();
    float& d();
    float& e();
    float& f();

    QVector3D abc()const;
    QVector3D def() const;


    const Quadric2D& operator=(const Quadric2D& q);

};

#endif // QUADRIC2D_H
