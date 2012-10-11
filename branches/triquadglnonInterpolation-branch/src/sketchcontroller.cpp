#include "sketchcontroller.h"
#include <QtOpenGL>

SketchController::SketchController(QObject *parent) :
    QObject(parent)
{
}


void SketchController::mouseRigthMove(const QPoint &curr, const QVector4D &proj)
{
    curve2.append(curr);
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
    return curve2;
}
