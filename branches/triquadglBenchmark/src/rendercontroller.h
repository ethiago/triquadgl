#ifndef RENDERCONTROLLER_H
#define RENDERCONTROLLER_H

#include <QObject>
#include <QPoint>
#include <QGLContext>
#include <QAction>
#include <QMap>
#include <QList>
#include <QPair>

class MainWindow;
class Object3D;
class TriQuadMesh;

class GLDisplay;
class SketchController;

class RenderController : public QObject
{
    Q_OBJECT

    GLDisplay *display;
    TriQuadMesh *triquad;
    SketchController *skC;

public:
    explicit RenderController(MainWindow *mainWindow,
                              QObject *parent = 0);
    ~RenderController();
    QGLWidget* getGLContext(void);
    void updateGL(void);

private:
    QString saveImage();

public slots:
    void drawModel(void);
    void mouseRigthMove(QPoint ini, QPoint curr);
    void mouseRigthFinish(QPoint ini, QPoint curr);
    void mouseLeftMove(QPoint ini, QPoint curr);
    void mouseLefthFinish(QPoint ini, QPoint curr);
    void mouseCancel();
    void saveResultAsImage(const QString& fn = QString());
    void viewMesh(bool);
    void benchmark();
};

#endif // RENDERCONTROLLER_H
