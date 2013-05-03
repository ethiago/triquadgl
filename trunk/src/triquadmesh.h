#ifndef TRIQUADMESH_H
#define TRIQUADMESH_H

#include "QGLShaderProgram"
#include "Object3D.h"
#include "QVector"
#include "quadric2d.h"
#include "compacthalfedge.h"

typedef struct
{
    int idx[3];
}NO;

class TriQuadMesh : public Object3D
{
    //QVector<QVector4D> inPoints;
    //QVector<QVector2D> vertices;
    //QVector<Quadric2D> quadrics;
    CompactHalfEdge triquads;

    QGLShaderProgram program;
    int locationABC;
    int locationDEF;

    bool showMesh;
    bool showSketch;
    bool origin;

    int vwp[4];
    QMatrix4x4 mvpi;

    int idxMaisProximo;
    QVector2D maisProximo;

    QVector<QVector2D> dp;
    QVector<QVector2D> dpL;
    QVector<QVector2D> dpU;

public:
    explicit TriQuadMesh(const QVector3D& center = QVector3D(),
                    QObject *parent = 0);
    explicit TriQuadMesh(const TriQuadMesh& tt);
    ~TriQuadMesh();

    virtual Object3D* copy() const;

    bool isProgramLinked();
    void viewMesh(bool v);
    void viewSketch(bool v);

    void changeOrigin(bool v);

    QVector4D unproject(const QPoint&);

    void move(const QPoint& ini, const QPoint& curr);
    void finish();
    void cancel();
    void fitting(const QVector<QPoint> &);
    void fittingG(QVector<QVector4D> &in);
    void globalFitting_2layers(QVector<QVector4D> & in);
    void globalFitting_2layers_freef(QVector<QVector4D> & in);
    void globalFittingG_3layers_freef(QVector<QVector4D> & in);
    void globalFitting_3layers(QVector<QVector4D> & in);
    void globalFitting_5layers(QVector<QVector4D> & in);
    void fittingGG(QVector<QVector4D> &in);
    void globalFittingWithNormals(QVector<QVector4D> & in);
    void addVertex(const QVector4D& newVertex);
    void joinVerticesAt(const QVector4D& controlPoint);

    void clearDrawPoints();

private:
    void drawOrigin();
    void drawPoints(const QVector<QVector2D>& ps = QVector<QVector2D>());
    void drawPoints1(const QVector<QVector2D>& ps = QVector<QVector2D>());
    void drawPoints2(const QVector<QVector2D>& ps = QVector<QVector2D>());
    virtual void drawGeometry(void);
    virtual void beforeTransformations(void);
    virtual void afterTransformations(void);
    QMatrix4x4 glGetMatrix(GLenum fetchType);
    void buildObject();
    QMatrix4x4 buildInv(int triangleId);
    QVector4D buscaNo(const QVector4D& p, int *idx);
    int configPoints(QVector<QVector4D>& pontos, QVector<QVector3D>& bary, QVector<int>& idx );
    QVector3D bary(const QVector4D& p, int *idx);
};

QDebug operator<< (QDebug d, const Quadric2D &model);


#endif // TRIQUADMESH_H
