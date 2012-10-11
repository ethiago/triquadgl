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
    connect(ui->actionSave_Mesh, SIGNAL(triggered()), this, SIGNAL(saveMesh()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setGLDisplay(GLDisplay *display)
{
    ui->verticalLayout->addWidget(display);
}
