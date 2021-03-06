#ifndef SKETCHCONTROLLER_H
#define SKETCHCONTROLLER_H

#include <QObject>
#include <QVector>
#include <QVector4D>

#define SKETCHFILEEXTENSION "skc"

class TriQuadMesh;

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

    bool loadSketchScreenCoordinates(const QString& fileName);
    bool loadSketch(const QString& fileName, TriQuadMesh *triquad);
    bool saveSketch(const QString& fileName, TriQuadMesh* triquad);


};



#endif // SKETCHCONTROLLER_H
