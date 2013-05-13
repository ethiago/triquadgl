#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
    class MainWindow;
}

class GLDisplay;

class MainWindow : public QMainWindow
{
    Q_OBJECT

signals:
    void saveResultAsImage();
    void viewMesh(bool);
    void loadMesh();
    void metodoMudou(int);
    void cheBuilderMudou(int);
    void loadSketch();
    void viewSketch(bool);
    void saveMesh();
    void saveSketch();
    void viewScalarField(bool);
    void clearMesh();
    void linearFilter(bool);
    void configsUpdated();

public:

    typedef struct _GRIDOPTIONS
    {
        float xm, xM, ym, yM;
        int xN, yN;
    }GRIDOPTIONS;

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void setGLDisplay(GLDisplay *);
    int metodoSelecionado();
    int cheBuilderSelecionado();
    void addMetodo(const QString&);
    void addCHEBuilder(const QString& label);
    GRIDOPTIONS getGridOptions();
    float getKDistance();
    bool includeVertices();

private:
    Ui::MainWindow *ui;

};

#endif // MAINWINDOW_H
