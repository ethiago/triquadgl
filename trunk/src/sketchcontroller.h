#ifndef SKETCHCONTROLLER_H
#define SKETCHCONTROLLER_H

#include <QObject>
#include <QVector>
#include <QVector4D>

#define SKETCHFILEEXTENSION "skc"

class SketchController : public QObject
{
    Q_OBJECT

    QVector<QPoint> last;
    QVector<QPoint> curve;
    QVector<QVector4D> feedBackPoints;

public:
    explicit SketchController(QObject *parent = 0);
    void draw();
    void mouseRigthMove(const QPoint& curr, const QVector4D& proj);
    void mouseRigthFinish();
    void cancel();

    QVector<QPoint> getPoints(void);
    QVector<QPointF> getPointsLinearFilter(void);

    bool loadSketch(const QString& fileName);
    bool saveSketch(const QString& fileName);


};



#endif // SKETCHCONTROLLER_H
