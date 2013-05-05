#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "GLDisplay.h"
#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->actionSaveImageAs, SIGNAL(triggered()), this, SIGNAL(saveResultAsImage()));
    connect(ui->actionView_Mesh, SIGNAL(toggled(bool)), this, SIGNAL(viewMesh(bool)));
    connect(ui->actionOpen_Mesh, SIGNAL(triggered()), this, SIGNAL(loadMesh()));
    connect(ui->cmb_Metodo, SIGNAL(activated(int)), SIGNAL(metodoMudou(int)));
    connect(ui->actionLoad_Sketch, SIGNAL(triggered()), this, SIGNAL(loadSketch()) );
    connect(ui->actionView_Sketch, SIGNAL(toggled(bool)), this, SIGNAL(viewSketch(bool)));
    connect(ui->actionSave_Mesh, SIGNAL(triggered()), this, SIGNAL(saveMesh()));
    connect(ui->actionSave_Sketch, SIGNAL(triggered()), this, SIGNAL(saveSketch()));
    connect(ui->actionView_Scalar_Field, SIGNAL(toggled(bool)),this, SIGNAL(viewScalarField(bool)));
    connect(ui->actionClear_Mesh, SIGNAL(triggered()), this, SIGNAL(clearMesh()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setGLDisplay(GLDisplay *display)
{
    ui->verticalLayout->addWidget(display);
}

int MainWindow::metodoSelecionado()
{
    return ui->cmb_Metodo->currentIndex();
}

void MainWindow::addMetodo(const QString& label)
{
    ui->cmb_Metodo->addItem(label);
}
