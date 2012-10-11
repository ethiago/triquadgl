#include <QtGui/QApplication>
#include <QDebug>
#include "mainwindow.h"
#include "rendercontroller.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow mw;
    RenderController rc(&mw);

//    const char *extensionString = reinterpret_cast<const char *>(glGetString(GL_EXTENSIONS));
//    const char *renderer = reinterpret_cast<const char *>(glGetString(GL_RENDERER));

//    if(extensionString == 0)
//        return 0;

//    QString rend(renderer);
//    qDebug() << rend;
//    QStringList extension = QString( extensionString).split(QChar(' '));
//    qDebug() << extension;

    return a.exec();
}
