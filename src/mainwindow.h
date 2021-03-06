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
    void loadSketchSC();
    void viewSketch(bool);
    void saveMesh();
    void saveSketch();
    void viewScalarField(bool);
    void clearMesh();
    void linearFilter(bool);
    void configsUpdated();
    void fittingMeasure();
    void meshTranslation(bool);
    void viewGradField(bool);
    void showTriQuad(bool);
    void lengthForVis(bool);
    void makeSmooth(void);
    void isoform(bool);
    void imageOpened(const QImage& );
    void localTriQuad(bool);

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

    bool isMeshView();
    bool isFieldView();
    bool isSketchView();

    void setStatusText(const QString& text);

public slots:
    void openImage();

private:
    Ui::MainWindow *ui;

};

#endif // MAINWINDOW_H
