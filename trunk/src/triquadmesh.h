#ifndef TRIQUADMESH_H
#define TRIQUADMESH_H

#include "QGLShaderProgram"
#include "Object3D.h"
#include "QVector"
#include "quadric2d.h"
#include "compacthalfedge.h"
#include "chebuilder.h"

#define MESHFILEEXTENSION "msh"

typedef struct
{
    int idx[3];
}NO;

class TriQuadMesh : public Object3D
{
    CompactHalfEdge che;

    QGLShaderProgram program[2];
    int locationABC[2];
    int locationDEF[2];
    int locationScalar[2];
    int locationMenor[2];
    int locationMaior[2];
    int locationTexture[2];
    int locationQxABC[2];
    int locationQxDEF[2];
    int locationQyABC[2];
    int locationQyDEF[2];
    int activeProgram;

    bool showMesh;
    bool showSketch;
    bool origin;
    bool showScalarField;
    bool showInputLine;
    bool showTriQuad;
    bool meshTranslation;

    int vwp[4];
    QMatrix4x4 mvpi;

    int idxMaisProximo;
    QVector2D maisProximo;

    QVector<QVector2D> dp;
    QVector<QVector2D> dpL;
    QVector<QVector2D> dpU;

    int m_glTextureName;

public:
    explicit TriQuadMesh(int texName, const QVector3D& center = QVector3D(),
                    QObject *parent = 0);
    explicit TriQuadMesh(const TriQuadMesh& tt);
    ~TriQuadMesh();

    virtual Object3D* copy() const;

    void clear();
    bool isEmpty();
    bool isProgramLinked();
    void viewMesh(bool v);
    void viewSketch(bool v);
    void viewGrad(bool v);

    void changeOrigin(bool v);

    void setMeshTranslation(bool);

    QVector4D unproject(const QPoint&);
    QVector4D unproject(const QPointF&);
    QVector<QVector4D> unproject(const QVector<QPoint>&);
    QVector<QVector4D> unproject(const QVector<QPointF>&);

    void move(const QPoint& ini, const QPoint& curr);
    void finish();
    void cancel();

    void globalFitting_2layers(QVector<QVector4D> in, float k);
    void globalFitting_2layers_freef(QVector<QVector4D>  in, float k);
    void globalFittingG_3layers_freef(QVector<QVector4D>  in, float k, bool includeVertices);
    void globalFittingG_3layers_freef_withGrad(QVector<QVector4D> pontos, float kDistance, bool includeVertices);
    void globalFittingG_3layers_freef_kDistance(QVector<QVector4D> pontos, float kDistance, bool includeVertices);
    void globalFittingG_1layers_freef(QVector<QVector4D> pontos);
    void globalFitting_3layers(QVector<QVector4D>  in, float k, bool includeVertices);
    void globalFitting_5layers(QVector<QVector4D>  in, float k);
    void globalFitting_1layer(QVector<QVector4D>  in, bool includeVertices);

    void addVertex(const QVector4D& newVertex);
    void joinVerticesAt(const QVector4D& controlPoint);
    void deleteTriangleWith(const Vertex& v);
    void viewScalarField(bool v);
    bool buildMesh(CHEBuilder* builder);

    void clearDrawPoints();
    bool loadMesh(const QString &filename);
    bool saveMesh(const QString &filename)const;

    void configureRenderInputLine();
    void configureRenderTriQuad();
    void reconfigure(bool,bool,bool);

    int textureName()const;

    QVector<QVector2D> pointsOnEdge(const HalfEdge &h, int c = 3)const;



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
    void calcDistances(const QVector<QVector4D> &pontos, QVector<float>* distances, QVector<int> *idx);

public slots:
    void resizeWindow(const QSize&);

};

QDebug operator<< (QDebug d, const Quadric2D &model);


#endif // TRIQUADMESH_H
