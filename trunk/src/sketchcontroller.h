#ifndef SKETCHCONTROLLER_H
#define SKETCHCONTROLLER_H

#include <QObject>
#include <QVector>
#include <QVector4D>

class SketchController : public QObject
{
    Q_OBJECT

    QVector<QVector4D> curve;

public:
    explicit SketchController(QObject *parent = 0);
    void draw();
    void mouseRigthMove(const QVector4D& curr);
    void mouseRigthFinish();
    void cancel();
    const QVector<QVector4D>& getPoints(void);

private:
    void processaCurva();

};



#endif // SKETCHCONTROLLER_H
