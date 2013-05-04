#include "sketchcontroller.h"
#include <QtOpenGL>

SketchController::SketchController(QObject *parent) :
    QObject(parent)
{
}


void SketchController::mouseRigthMove(const QPoint &curr, const QVector4D &proj)
{
    curve.append(curr);
    feedBackPoints.append(proj);
}

void SketchController::mouseRigthFinish()
{
    last = curve;
    curve.clear();
    feedBackPoints.clear();
}

void SketchController::cancel()
{
    curve.clear();
    feedBackPoints.clear();
}

void SketchController::draw()
{
    if(feedBackPoints.size() == 0)
        return;

    glBegin(GL_POINTS);
    {
        for(int i = 0; i < feedBackPoints.size(); ++i)
        {
            glVertex4fv(reinterpret_cast<const GLfloat *>(&feedBackPoints[i]));
        }
    }glEnd();
}

QVector<QPoint> SketchController::getPoints(void)
{
    return last;
}

bool SketchController::loadSketch(const QString& fn )
{
    QFile f(fn);

    if(!f.open(QIODevice::ReadOnly))
        return false;

    QTextStream s(f.readAll());
    f.close();

    last.clear();

    int qtd;
    s >> qtd;

    for(int i = 0; i < qtd; ++i)
    {
        int x,y;
        s >> x >> y;
        last.append(QPoint(x,y));
    }
    return true;
}

bool SketchController::saveSketch(const QString& fileName )
{
    QFile f(fileName);

    if(!f.open(QIODevice::WriteOnly))
        return false;

    QTextStream s(&f);

    s << last.size() << "\n";
    for(int i = 0; i < last.size(); ++i)
    {
        s << last[i].x() << " " << last[i].y() << "\n";
    }

    f.close();
    return true;
}
