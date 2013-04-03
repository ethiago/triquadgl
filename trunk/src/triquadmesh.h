#ifndef TRIQUADMESH_H
#define TRIQUADMESH_H

#include "QGLShaderProgram"
#include "Object3D.h"
#include "QVector"

#define CIRCLE   TriQuadMesh::makeQuadric(1.0, 1.0, 0.0, 0.0, 0.0,-0.3)
#define CIRCLE2  TriQuadMesh::makeQuadric(1.0, 1.0, 0.0,-2.0,-4.0, 4.0)
#define PARABOLA TriQuadMesh::makeQuadric(2.0, 0.0, 0.0, 0.0, 1.0, 0.0)
#define ZERO     TriQuadMesh::makeQuadric(0.0, 0.0, 0.0, 0.0, 0.0, 1.0)

typedef struct
{
    QVector3D a_b_c;
    QVector3D d_e_f;
}Quadric;

typedef struct
{
    int idx[3];
}NO;

class TriQuadMesh : public Object3D
{
    //QVector<QVector4D> inPoints;
    QVector<QVector2D> vertices;
    QVector<Quadric> quadrics;
    QVector<NO> triquads;

    QGLShaderProgram program;
    int locationABC;
    int locationDEF;

    bool showMesh;
    bool origin;

    int vwp[4];
    QMatrix4x4 mvpi;

    int idxMaisProximo;
    QVector2D maisProximo;

public:
    explicit TriQuadMesh(const QVector3D& center = QVector3D(),
                    QObject *parent = 0);
    explicit TriQuadMesh(const TriQuadMesh& tt);
    ~TriQuadMesh();

    virtual Object3D* copy() const;

    void setQuadric(int idx, const Quadric& c);

    bool isProgramLinked();
    void viewMesh(bool v);

    static Quadric makeQuadric(float x2, float y2, float xy, float x, float y, float c);

    void changeOrigin(bool v);

    QVector4D unproject(const QPoint&);
    int busca(const QVector4D&);

    void move(const QPoint& ini, const QPoint& curr);
    void finish();
    void cancel();
    void fitting(const QVector<QPoint> &);
    void fittingG(QVector<QVector4D> &in);
    void fittingG2(QVector<QVector4D> & in);
    void fittingGG(QVector<QVector4D> &in);
    void globalFittingWithNormals(QVector<QVector4D> & in);

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
    QMatrix4x4 buildInv(NO&);
    QVector4D buscaNo(const QVector4D& p, int *idx);
    int configPoints(QVector<QVector4D>& pontos, QVector<QVector3D>& bary, QVector<int>& idx );
    QVector3D bary(const QVector4D& p, int *idx);
};

QDebug operator<< (QDebug d, const Quadric &model);


#endif // TRIQUADMESH_H
