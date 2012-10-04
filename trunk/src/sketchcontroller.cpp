#include "sketchcontroller.h"
#include <QtOpenGL>

SketchController::SketchController(QObject *parent) :
    QObject(parent)
{
}


void SketchController::mouseRigthMove(const QVector4D &curr)
{
    curve.append(curr);
}

void SketchController::mouseRigthFinish()
{
    processaCurva();
}

void SketchController::cancel()
{
    curve.clear();
}

void SketchController::processaCurva()
{
}

void SketchController::draw()
{
    if(curve.size() == 0)
        return;

    glBegin(GL_POINTS);
    {
        for(int i = 0; i < curve.size(); ++i)
        {
            glVertex4fv(reinterpret_cast<const GLfloat *>(&(curve[i])));
        }
    }glEnd();
}

const QVector<QVector4D> &SketchController::getPoints(void)
{
    return curve;
}
