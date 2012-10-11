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

QVector<QVector2D> SketchController::loadSketch(const QString& fn )
{
    QVector<QVector2D> ret;
    QFile f(fn);

    if(!f.open(QIODevice::ReadOnly))
        return ret;

    QTextStream s(f.readAll());
    f.close();

    int qtd;
    s >> qtd;

    for(int i = 0; i < qtd; ++i)
    {
        float x,y;
        s >> x >> y;
        ret.append(QVector2D(x,y));
    }
    return ret;
}
