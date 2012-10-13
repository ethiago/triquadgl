#ifndef CUBIC_H
#define CUBIC_H

#include <QVector3D>
#include <QGLShaderProgram>

#include "Object3D.h"

class Cubic : public Object3D
{
    QVector2D m_v[3];
    float consts[10];
    QGLShaderProgram program;
    GLint constsLocation;

public:
    Cubic();

    QVector2D& operator [](int i);

    void setConsts(float x3, float x2y, float xy2, float y3, float x2, float xy, float y2, float x, float y, float c );

    void bulild();

    float * getConsts();

private:
    virtual void beforeTransformations(void);
    virtual void afterTransformations(void);
    virtual void drawGeometry(void);

};

#endif // CUBIC_H
