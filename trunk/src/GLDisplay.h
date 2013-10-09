#ifndef GLDISPLAY_H
#define GLDISPLAY_H

#include <QGLWidget>
#include <QMouseEvent>
#include <QMoveEvent>
#include <QPoint>


#define NULLPOINT   QPoint(-1,-1)
#define EPSILON     0.001

class GLDisplay : public QGLWidget
{
    Q_OBJECT

    QPoint rigthPressedPoint;
    QPoint leftPressedPoint;
    float zoom;


signals:
    void drawModel(void);
    void mouseRigthMove(QPoint ini, QPoint curr);
    void mouseRigthFinish(QPoint ini, QPoint curr);
    void mouseLeftMove(QPoint ini, QPoint curr);
    void mouseLefthFinish(QPoint ini, QPoint curr);
    void mouseCancel();
    void mouseDoubleClickLeft(QPoint);
    void resizeWindow(const QSize&);

public:
    GLDisplay(QWidget *parent = 0);
    ~GLDisplay();

    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

private:
    void mousePressEvent ( QMouseEvent * event );
    void mouseReleaseEvent ( QMouseEvent * event );
    void mouseMoveEvent(QMouseEvent *);
    void wheelEvent(QWheelEvent * event);
    void mouseDoubleClickEvent(QMouseEvent *);
    float xDist(float aspect);
    float yDist(float aspect);

};


#endif // GLDISPLAY_H
