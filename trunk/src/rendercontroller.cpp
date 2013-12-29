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
#include "chebuilderdefault.h"
#include "chebuilderregulargrid.h"
#include "chebuilderquadtreefrompointcloud.h"
#include "chebuilderregulargridfrompointcloud.h"
#include "chebuilderequilateralmesh.h"
#include "fastmarching.h"

RenderController::RenderController(MainWindow *mainWindow,
                                   QObject *parent):
    QObject(parent)
{
    mw = mainWindow;
    skC = new SketchController();
    this->display = new GLDisplay();
    mainWindow->setGLDisplay(display);
    metodo = 0;
    cheBuilder = 0;
    m_linearFilter = false;
    configComboMetodo();
    configComboCHEBuilder();
    {  // esta ordem deve ser mantida
        display->updateGL();

        int textureName = display->bindTexture(QImage(":/texpontos"), GL_TEXTURE_2D);

        triquad = new TriQuadMesh(textureName);

        connect(display, SIGNAL(drawModel()),
                this, SLOT(drawModel()));
        connect(display, SIGNAL(resizeWindow(QSize)),
                triquad, SLOT(resizeWindow(QSize)));
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

    connect(mainWindow,SIGNAL(cheBuilderMudou(int)),
            this, SLOT(cheBuilderMudou(int)));

    connect(mainWindow, SIGNAL(viewScalarField(bool)),
            this, SLOT(viewScalarField(bool)));

    connect(mainWindow, SIGNAL(clearMesh()),
            this, SLOT(clearMesh()));

    connect(mainWindow, SIGNAL(linearFilter(bool)),
            this, SLOT(linearFilter(bool)) );

    connect(mainWindow, SIGNAL(configsUpdated()),
            this, SLOT(configsUpdated()) );

    connect(mainWindow, SIGNAL(fittingMeasure()),
            this, SLOT(fittingMeasure()) );

    connect(mainWindow, SIGNAL(meshTranslation(bool)),
            this, SLOT(meshTranslation(bool)) );

    connect(mainWindow, SIGNAL(viewGradField(bool)),
            this, SLOT(viewGradField(bool)) );

    connect(mainWindow, SIGNAL(showTriQuad(bool)),
            this, SLOT(showTriQuad(bool)) );

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
    ini.rx();
    skC->mouseRigthMove(curr, triquad->unproject(curr));
    display->updateGL();
}

void RenderController::mouseRigthFinish(QPoint ini, QPoint curr)
{
    if(ini != curr)
    {
        skC->mouseRigthFinish();
        if(m_linearFilter)
            ultimaLista = triquad->unproject(skC->getPointsLinearFilter());
        else
            ultimaLista = triquad->unproject(skC->getPoints());
        exec();
    }else
    {
        triquad->deleteTriangleWith(triquad->unproject(curr).toVector2D());
        display->updateGL();
    }
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

    skC->loadSketch(filename, triquad);

    if(m_linearFilter)
        ultimaLista = triquad->unproject(skC->getPointsLinearFilter());
    else
        ultimaLista = triquad->unproject(skC->getPoints());
    exec();
}

void RenderController::metodoMudou(int m )
{
    metodo = m;
    exec();
}

void RenderController::cheBuilderMudou(int v )
{
    cheBuilder = v;
}

void RenderController::clearMesh()
{
    triquad->clear();
    display->updateGL();
}

void RenderController::configComboMetodo()
{
    mw->addMetodo("1 Camada");
    mw->addMetodo("3 Camadas");
    mw->addMetodo("2 Camadas");
    mw->addMetodo("5 Camadas");
    mw->addMetodo("2 Camadas - f livre");
    mw->addMetodo("3 Camadas - f livre");
    mw->addMetodo("3 Camadas - f livre (K Distance)");
    mw->addMetodo("1 Camada  - f livre (K Distance)");
    mw->addMetodo("3 Camadas - f livre - grad");
    mw->addMetodo("Q");
    metodo = mw->metodoSelecionado();
}

void RenderController::configComboCHEBuilder()
{
    mw->addCHEBuilder("Default");
    mw->addCHEBuilder("Regular Grid");
    mw->addCHEBuilder("Regular Grid From Point Cloud");
    mw->addCHEBuilder("QuadTree From Point Cloud");
    mw->addCHEBuilder("Equilateral Mesh From Point Cloud");
    metodo = mw->metodoSelecionado();
}

void RenderController::exec()
{
    if(triquad->isEmpty())
    {
        MainWindow::GRIDOPTIONS opt = mw->getGridOptions();
        CHEBuilder *builder = NULL;
        switch (cheBuilder)
        {
        case 0:
            builder = new CHEBuilderDefault();
            break;
        case 1:
            builder = new CHEBuilderRegularGrid(opt.xm, opt.xM, opt.ym, opt.yM, opt.xN, opt.yN);
            break;
        case 2:
            builder = new CHEBuilderRegularGridFromPointCloud(ultimaLista, opt.xN, opt.yN);
            break;
        case 3:
            builder = new CHEBuilderOctreeFromPointCloud(ultimaLista);
            break;
        case 4:
            builder = new CHEBuilderEquilateralMesh(ultimaLista, opt.xN);
            break;
        }
        bool built = triquad->buildMesh(builder);
        if(builder)
            delete builder;
        if(!built)
            return;
    }

    float k = mw->getKDistance();
    bool includeVertices = mw->includeVertices();

    switch (metodo)
    {
    case 0:
        triquad->globalFitting_1layer(ultimaLista, includeVertices);
        break;
    case 1:
        triquad->globalFitting_3layers(ultimaLista, k, includeVertices);
        break;
    case 2:
        triquad->globalFitting_2layers(ultimaLista, k);
        break;
    case 3:
        triquad->globalFitting_5layers(ultimaLista, k);
        break;
    case 4:
        triquad->globalFitting_2layers_freef(ultimaLista, k);
        break;
    case 5:
        triquad->globalFittingG_3layers_freef(ultimaLista, k, includeVertices);
        break;
    case 6:
        triquad->globalFittingG_3layers_freef_kDistance(ultimaLista, k, includeVertices);
        break;
    case 7:
        triquad->globalFittingG_1layers_freef(ultimaLista);
        break;
    case 8:
        triquad->globalFittingG_3layers_freef_withGrad(ultimaLista, k, includeVertices);
        break;
    case 9:
        triquad->fitting_quadrica(ultimaLista);
        break;

    }
    //fittingMeasure();
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

    skC->saveSketch(filename, triquad);
}

void RenderController::viewScalarField(bool v)
{
    triquad->viewScalarField(v);
    display->updateGL();
}

void RenderController::linearFilter(bool v)
{
    m_linearFilter = v;
    if(m_linearFilter)
        ultimaLista = triquad->unproject(skC->getPointsLinearFilter());
    else
        ultimaLista = triquad->unproject(skC->getPoints());
    exec();
}

void RenderController::configsUpdated()
{
    exec();
}

void RenderController::fittingMeasure()
{
    triquad->configureRenderInputLine();
    display->updateGL();
    display->updateGL();
    QImage input = display->grabFrameBuffer();

    triquad->configureRenderTriQuad();
    display->updateGL();
    display->updateGL();
    QImage result = display->grabFrameBuffer();

    triquad->reconfigure(mw->isMeshView(), mw->isSketchView(), mw->isFieldView());
    display->updateGL();
    display->updateGL();

    input.save("srcInput.png");
    FastMarching fm(input);
    fm.run();
    fm.getImage().save("input.png");

    result.save("srcTriquad.png");
    FastMarching fm2(result);
    fm2.run();
    fm2.getImage().save("triquad.png");

    float howFitted   = fm.distanceTo(fm2);
    float moreThenFit = fm2.distanceTo(fm);

    mw->setStatusText(QString::number(howFitted) + " X " + QString::number(moreThenFit));
    mw->update();
}

void RenderController::meshTranslation(bool v)
{
    triquad->setMeshTranslation(v);
}

void RenderController::viewGradField(bool v)
{
    triquad->viewGrad(v);
    display->updateGL();
}

void RenderController::showTriQuad(bool v)
{
    triquad->viewTriQuad(v);
    display->updateGL();
}
