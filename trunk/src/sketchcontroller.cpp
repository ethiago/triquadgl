#include "sketchcontroller.h"
#include <QtOpenGL>
#include "vectorn.hpp"

SketchController::SketchController(QObject *parent) :
    QObject(parent)
{
}


void SketchController::mouseRigthMove(const QPoint &curr, const QVector4D &proj)
{
    PointN<float,2> p;
    p[0] = curr.x();
    p[1] = curr.y();
    curve2.add(p);
    curve.append(proj);
}

void SketchController::mouseRigthFinish()
{
    processaCurva();
}

void SketchController::cancel()
{
    curve2.clear();
    curve.clear();
}

void SketchController::processaCurva()
{
   // curve2.lineFilter();
}

void SketchController::draw()
{
    if(curve.size() == 0)
        return;

    glBegin(GL_POINTS);
    {
        for(int i = 0; i < curve.size(); ++i)
        {
            glVertex4fv(reinterpret_cast<const GLfloat *>(&curve[i]));
        }
    }glEnd();
}

QVector<QPoint> SketchController::getPointsLinearFilter(void)
{
    QVector<QPoint> resp;
    for(int i = 0; i < curve2.size(); ++i)
    {
        resp.push_back(QPoint(curve2[i][0],curve2[i][1]));
    }
    return resp;
}
