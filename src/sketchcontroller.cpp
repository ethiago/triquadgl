#include "sketchcontroller.h"
#include <QtOpenGL>
#include "curveN/curven.hpp"
#include "triquadmesh.h"

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

QVector<QPointF> SketchController::getPointsLinearFilter(void)
{
    PolygonalCurve <float , 2> curveFilter;
    for(int i = 0; i < last.size(); ++i)
    {
        curveFilter.add(PointN<float,2>(last[i].x(), last[i].y()));
    }
    curveFilter.lineFilter(0.5,4);

    QVector<QPointF> ret;
    for(unsigned int i = 0; i < curveFilter.size(); ++i)
    {
        ret.append(QPointF(curveFilter[i][0], curveFilter[i][1]));
    }
    return ret;

}

bool SketchController::loadSketch(const QString& fn , TriQuadMesh* triquad)
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
        float x,y;
        s >> x >> y;

        QPoint p = triquad->project(QVector4D(x,y,0.0,1.0));
        qDebug() << p;

        last.append(p);
    }
    return true;
}

bool SketchController::loadSketchScreenCoordinates(const QString& fn)
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

bool SketchController::saveSketch(const QString& fileName, TriQuadMesh* triquad )
{
    QFile f(fileName);

    if(!f.open(QIODevice::WriteOnly))
        return false;

    QTextStream s(&f);

    s << last.size() << "\n";
    for(int i = 0; i < last.size(); ++i)
    {
        QVector4D v = triquad->unproject(last[i]);
        s << v.x() << " " << v.y() << "\n";
    }

    f.close();
    return true;
}
