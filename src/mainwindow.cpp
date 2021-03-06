#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "GLDisplay.h"
#include <QDebug>
#include <QFileDialog>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->actionSaveImageAs, SIGNAL(triggered()), this, SIGNAL(saveResultAsImage()));
    connect(ui->actionView_Mesh, SIGNAL(toggled(bool)), this, SIGNAL(viewMesh(bool)));
    connect(ui->actionOpen_Mesh, SIGNAL(triggered()), this, SIGNAL(loadMesh()));
    connect(ui->cmb_Metodo, SIGNAL(activated(int)), SIGNAL(metodoMudou(int)));
    connect(ui->cmb_cheBuilder, SIGNAL(activated(int)), this, SIGNAL(cheBuilderMudou(int)));
    connect(ui->actionLoad_Sketch, SIGNAL(triggered()), this, SIGNAL(loadSketch()) );
    connect(ui->actionLoad_From_Screen_Coordinate, SIGNAL(triggered()), this, SIGNAL(loadSketchSC()) );
    connect(ui->actionView_Sketch, SIGNAL(toggled(bool)), this, SIGNAL(viewSketch(bool)));
    connect(ui->actionSave_Mesh, SIGNAL(triggered()), this, SIGNAL(saveMesh()));
    connect(ui->actionSave_Sketch, SIGNAL(triggered()), this, SIGNAL(saveSketch()));
    connect(ui->actionView_Scalar_Field, SIGNAL(toggled(bool)),this, SIGNAL(viewScalarField(bool)));
    connect(ui->actionClear_Mesh, SIGNAL(triggered()), this, SIGNAL(clearMesh()));
    connect(ui->actionLinear_Filter, SIGNAL(toggled(bool)), this, SIGNAL(linearFilter(bool)));
    connect(ui->kDistance, SIGNAL(editingFinished()), this, SIGNAL(configsUpdated()) );
    connect(ui->includeVertices, SIGNAL(clicked()), this, SIGNAL(configsUpdated()) );
    connect(ui->actionFitting_Measure, SIGNAL(triggered()), this, SIGNAL(fittingMeasure()) );
    connect(ui->meshTrans, SIGNAL(toggled(bool)), this, SIGNAL(meshTranslation(bool)) );
    connect(ui->actionView_Grad_Field, SIGNAL(toggled(bool)), this, SIGNAL(viewGradField(bool)) );
    connect(ui->actionShow_TriQuad, SIGNAL(toggled(bool)), this, SIGNAL(showTriQuad(bool)) );
    connect(ui->lengthVis, SIGNAL(toggled(bool)), this, SIGNAL(lengthForVis(bool)) );
    connect(ui->actionMake_Smooth, SIGNAL(triggered()), this, SIGNAL(makeSmooth()) );
    connect(ui->iso, SIGNAL(toggled(bool)), this, SIGNAL(isoform(bool)) );
    connect(ui->actionLoad_Image, SIGNAL(triggered()), this, SLOT(openImage()) );
    connect(ui->actionLocal_TriQuad, SIGNAL(toggled(bool)), this,  SIGNAL(localTriQuad(bool)) );

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setGLDisplay(GLDisplay *display)
{
    ui->renderLayout->addWidget(display);
}

int MainWindow::metodoSelecionado()
{
    return ui->cmb_Metodo->currentIndex();
}

int MainWindow::cheBuilderSelecionado()
{
    return ui->cmb_cheBuilder->currentIndex();
}

void MainWindow::addMetodo(const QString& label)
{
    ui->cmb_Metodo->addItem(label);
}

void MainWindow::addCHEBuilder(const QString& label)
{
    ui->cmb_cheBuilder->addItem(label);
}

MainWindow::GRIDOPTIONS MainWindow::getGridOptions()
{
    GRIDOPTIONS opt;
    opt.xm = ui->xMin->value();
    opt.xM = ui->xMax->value();
    opt.ym = ui->yMin->value();
    opt.yM = ui->yMax->value();
    opt.xN = ui->xBuckets->value();
    opt.yN = ui->yBuckets->value();
    return opt;
}

float MainWindow::getKDistance()
{
    return ui->kDistance->value();
}

bool MainWindow::includeVertices()
{
    return ui->includeVertices->isChecked();
}

bool MainWindow::isMeshView()
{
    return ui->actionView_Mesh->isChecked();
}

bool MainWindow::isFieldView()
{
    return ui->actionView_Scalar_Field->isChecked();
}

bool MainWindow::isSketchView()
{
    return ui->actionView_Sketch->isChecked();
}

void MainWindow::setStatusText(const QString& text)
{
    ui->statusBar->showMessage(text);
}

void MainWindow::openImage()
{
    QString fn = QFileDialog::getOpenFileName(this, "Open Image to fit", "..", tr("Images (*.png *.xpm *.jpg *.bmp *.gif)"));

    if(fn.isEmpty())
        return;

    QImage img(fn);
    if(img.isNull())
        return;

    emit imageOpened(img);
}
