#include <QDebug>
#include <QtOpenGL>
#include <QPoint>
#include <QFileDialog>

#include "rendercontroller.h"
#include "GLDisplay.h"
#include "mainwindow.h"
#include "Object3D.h"
#include "triquadmesh.h"
#include "sketchcontroller.h"

RenderController::RenderController(MainWindow *mainWindow,
                                   QObject *parent):
    QObject(parent)
{
    skC = new SketchController();
    this->display = new GLDisplay();
    mainWindow->setGLDisplay(display);
    metodo = 0;
    configCombo(mainWindow);
    {  // esta ordem deve ser mantida
        display->updateGL();

        triquad = new TriQuadMesh();

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

    connect(display, SIGNAL(mouseDoubleClickLeft(QPoint)),
            this, SLOT(mouseDoubleClickLeft(QPoint)));

    connect(mainWindow, SIGNAL(saveResultAsImage()),
            this, SLOT(saveResultAsImage()));

    connect(mainWindow, SIGNAL(viewMesh(bool)),
            this, SLOT(viewMesh(bool)));

    connect(mainWindow, SIGNAL(loadMesh()),
            this, SLOT(loadMesh()));

    connect(mainWindow, SIGNAL(loadSketch()),
            this, SLOT(loadSketch()));

    connect(mainWindow, SIGNAL(viewSketch(bool)),
            this, SLOT(viewSketch(bool)) );

    connect(mainWindow, SIGNAL(saveSketch()),
            this, SLOT(saveSketch()) );

    connect(mainWindow, SIGNAL(saveMesh()),
            this, SLOT(saveMesh()) );

    connect(mainWindow, SIGNAL(metodoMudou(int)),
            this, SLOT(metodoMudou(int)));

    connect(mainWindow, SIGNAL(viewScalarField(bool)),
            this, SLOT(viewScalarField(bool)));

    connect(mainWindow, SIGNAL(clearMesh()), this, SLOT(clearMesh()));

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
    triquad->draw();
    skC->draw();
}

void RenderController::mouseRigthMove(QPoint ini, QPoint curr)
{
    skC->mouseRigthMove(curr, triquad->unproject(curr));
    display->updateGL();
}

void RenderController::mouseRigthFinish(QPoint ini, QPoint curr)
{
    skC->mouseRigthFinish();
    ultimaLista = triquad->unproject(skC->getPoints());
    exec();
}

void RenderController::mouseLeftMove(QPoint ini, QPoint curr)
{
    triquad->move(ini,curr);
    display->updateGL();
}

void RenderController::mouseLefthFinish(QPoint ini, QPoint curr)
{
    triquad->finish();

    if(ini == curr)
    {
        temp = curr;
        timer.singleShot(1000, this, SLOT(timeout()) );
    }

    display->updateGL();
}

void RenderController::mouseDoubleClickLeft(QPoint p)
{
    if(temp != QPoint(-1,-1))
    {
        triquad->joinVerticesAt(triquad->unproject(p));
        temp = QPoint(-1,-1);
        timer.stop();
        display->updateGL();
    }
}

void RenderController::timeout()
{

    if(temp != QPoint(-1,-1))
    {
        triquad->addVertex(triquad->unproject(temp));
        temp = QPoint(-1,-1);
        display->updateGL();
    }
}

void RenderController::mouseCancel()
{
    skC->cancel();
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
    triquad->viewMesh(v);
    display->updateGL();
}

void RenderController::viewSketch(bool v)
{
    triquad->viewSketch(v);
    display->updateGL();
}

void RenderController::loadMesh()
{
    QString filename = QFileDialog::getOpenFileName(0, "Mesh File Loader", QString(), QString("*.")+MESHFILEEXTENSION );

    if(filename.isEmpty())
        return;

    triquad->loadMesh(filename);

    exec();
}

void RenderController::loadSketch()
{
    QString filename = QFileDialog::getOpenFileName(0, "Sketch File Loader", QString(), QString("*.")+SKETCHFILEEXTENSION );

    if(filename.isEmpty())
        return;

    skC->loadSketch(filename);

    ultimaLista = triquad->unproject(skC->getPoints());
    exec();
}

void RenderController::configCombo(MainWindow *mw)
{
    mw->addMetodo("3 Camadas");
    mw->addMetodo("2 Camadas");
    mw->addMetodo("5 Camadas");
    mw->addMetodo("2 Camadas - f livre");
    mw->addMetodo("3 Camadas - f livre");
    metodo = mw->metodoSelecionado();
}

void RenderController::metodoMudou(int m )
{
    metodo = m;
    qDebug() << metodo;
    exec();
}

void RenderController::clearMesh()
{
    triquad->clear();
    display->updateGL();
}

void RenderController::exec()
{
    if(triquad->isEmpty())
        triquad->buildMesh(ultimaLista);

    switch (metodo)
    {
    case 0:
        triquad->globalFitting_3layers(ultimaLista);
        break;
    case 1:
        triquad->globalFitting_2layers(ultimaLista);
        break;
    case 2:
        triquad->globalFitting_5layers(ultimaLista);
        break;
    case 3:
        triquad->globalFitting_2layers_freef(ultimaLista);
        break;
    case 4:
        triquad->globalFittingG_3layers_freef(ultimaLista);
        break;
    }
    display->updateGL();
}

void RenderController::saveMesh()
{
    QString filename = QFileDialog::getSaveFileName(0, "Save Mesh", QString(), QString(".")+MESHFILEEXTENSION );

    if(filename.isEmpty())
        return;

    triquad->saveMesh(filename);
}

void RenderController::saveSketch()
{
    QString filename = QFileDialog::getSaveFileName(0, "Save Sketch", QString(), QString(".")+SKETCHFILEEXTENSION );

    if(filename.isEmpty())
        return;

    skC->saveSketch(filename);
}

void RenderController::viewScalarField(bool v)
{
    triquad->viewScalarField(v);
    display->updateGL();
}