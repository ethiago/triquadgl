#ifndef SKETCHCONTROLLER_H
#define SKETCHCONTROLLER_H

#include <QObject>
#include <QVector>
#include <QVector4D>

class SketchController : public QObject
{
    Q_OBJECT

    QVector<QPoint> curve2;
    QVector<QVector4D> curve;

public:
    explicit SketchController(QObject *parent = 0);
    void draw();
    void mouseRigthMove(const QPoint& curr, const QVector4D& proj);
    void mouseRigthFinish();
    void cancel();
    QVector<QPoint> getPointsLinearFilter(void);

    static QVector<QVector2D> loadSketch(const QString& fn );

private:
    void processaCurva();

};



#endif // SKETCHCONTROLLER_H
