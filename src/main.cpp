//#include <QtWidgets>
#include <QtGui/QApplication>
#include <QDebug>
#include "mainwindow.h"
#include "rendercontroller.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow mw;
    RenderController rc(&mw);

    return a.exec();
}
