#include <QDebug>
#include <QtOpenGL>
#include <QPoint>
#include <QFileDialog>

#include "rendercontroller.h"
#include "GLDisplay.h"
#include "mainwindow.h"
#include "Object3D.h"
#include "triquadmesh.h"
#include "cubic.h"

RenderController::RenderController(MainWindow *mainWindow,
                                   QObject *parent):
    QObject(parent)
{
    this->display = new GLDisplay();
    mainWindow->setGLDisplay(display);
    {  // esta ordem deve ser mantida
        display->updateGL();

        triquad = new TriQuadMesh();
        cubic = new Cubic();

        cubic->setTraslation(QVector3D(-2,0,0));
        triquad->setTraslation(QVector3D(2,0,0));

        cubic->setConsts(1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0, 0.0);

        connect(display, SIGNAL(drawModel()),
                this, SLOT(drawModel()));
    }

    connect(display, SIGNAL(mouseLefthFinish(QPoint,QPoint)),
            this, SLOT(mouseLefthFinish(QPoint,QPoint)));

    connect(display, SIGNAL(mouseLeftMove(QPoint,QPoint)),
            this, SLOT(mouseLeftMove(QPoint,QPoint)));

    connect(display, SIGNAL(mouseRigthMove(QPoint,QPoint)),
            this, SLOT(mouseRigthMove(QPoint,QPoint)));

    connect(display, SIGNAL(mouseRigthFinish(QPoint,QPoint)),
            this, SLOT(mouseRigthFinish(QPoint,QPoint)));

    connect(display, SIGNAL(mouseCancel()),
            this, SLOT(mouseCancel()));

    connect(mainWindow, SIGNAL(saveResultAsImage()),
            this, SLOT(saveResultAsImage()));

    connect(mainWindow, SIGNAL(viewMesh(bool)),
            this, SLOT(viewMesh(bool)));

    connect(mainWindow, SIGNAL(buildTriQuad()),
            this, SLOT(buildTriQuad()));

    mainWindow->showMaximized();

}

RenderController::~RenderController()
{
    delete triquad;
    delete display;
}

void RenderController::updateGL(void)
{
    display->updateGL();
}


void RenderController::drawModel(void)
{
    triquad->getTriangle(0, (*cubic)[0], (*cubic)[1], (*cubic)[2]);
    triquad->draw();
    cubic->draw();
}

void RenderController::mouseRigthMove(QPoint ini, QPoint curr)
{
}

void RenderController::mouseRigthFinish(QPoint ini, QPoint curr)
{
}

void RenderController::mouseLeftMove(QPoint ini, QPoint curr)
{
    triquad->move(ini,curr);
    display->updateGL();
}

void RenderController::mouseLefthFinish(QPoint ini, QPoint curr)
{
    triquad->finish();
    display->updateGL();
}

void RenderController::mouseCancel()
{
    triquad->cancel();
    display->updateGL();
}

void RenderController::saveResultAsImage()
{
    QString filePath = saveImage();

    if(filePath.isEmpty())
        return;

    display->updateGL();
    display->grabFrameBuffer().save(filePath);
}

QGLWidget* RenderController::getGLContext(void)
{
    return display;
}

QString RenderController::saveImage()
{
    static QString dir = ".";
    QString filename = QFileDialog::getSaveFileName(0,
                                                    "Save Image",
                                                    dir,
                                                    "*.png");
    QString ext = ".png";

    if(filename == ext || filename.isEmpty())
        return QString();

    QFileInfo fi(filename);
    dir = fi.absolutePath();

    if(filename.right(ext.length()) != ext)
        filename += ext;

    return filename;
}

void RenderController::viewMesh(bool v)
{
    display->updateGL();
}

void RenderController::buildTriQuad()
{
    triquad->fromCubic(cubic->getConsts());
    display->updateGL();
}
