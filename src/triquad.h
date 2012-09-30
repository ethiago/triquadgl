#ifndef TRIQUAD_H
#define TRIQUAD_H

#include "QGLShaderProgram"
#include "Object3D.h"
#include "QVector"

#define CIRCLE   TriQuad::makeQuadric(1.0, 1.0, 0.0, 0.0, 0.0,-0.3)
#define CIRCLE2  TriQuad::makeQuadric(1.0, 1.0, 0.0,-2.0,-4.0, 4.0)
#define PARABOLA TriQuad::makeQuadric(2.0, 0.0, 0.0, 0.0, 1.0, 0.0)

typedef struct _Quadric
{
    QVector3D a_b_c;
    QVector3D d_e_f;
}Quadric;

typedef struct _NO
{
    int idx[3];
}NO;

class TriQuad : public Object3D
{
    QVector<QVector4D> vertices;
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
    QVector4D maisProximo;

public:
    explicit TriQuad(const QVector3D& center = QVector3D(),
                    QObject *parent = 0);
    explicit TriQuad(const TriQuad& tt);
    ~TriQuad();

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

private:
    void drawOrigin();
    virtual void drawGeometry(void);
    virtual void beforeTransformations(void);
    virtual void afterTransformations(void);
    QMatrix4x4 glGetMatrix(GLenum fetchType);
    void buildObject();
};

#endif // TRIQUAD_H
