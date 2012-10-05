#ifndef SKETCHCONTROLLER_H
#define SKETCHCONTROLLER_H

#include <QObject>
#include <QVector>
#include <QVector4D>
#include "curven.hpp"

typedef PolygonalCurve<float, 2> Curve;


class SketchController : public QObject
{
    Q_OBJECT

    Curve curve2;
    QVector<QVector4D> curve;

public:
    explicit SketchController(QObject *parent = 0);
    void draw();
    void mouseRigthMove(const QPoint& curr, const QVector4D& proj);
    void mouseRigthFinish();
    void cancel();
    QVector<QPoint> getPointsLinearFilter(void);

private:
    void processaCurva();

};



#endif // SKETCHCONTROLLER_H
