#include "triquadmesh.h"
#include <QGLShader>
#include <QVector3D>
#include <QtOpenGL>
#include "fitting.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "Curve.h"
#include "slgl3w.h"

#define SIGN(x) ((x)<0?-1.0:1.0)

TriQuadMesh::TriQuadMesh(int texName, const QVector3D& center, QObject *parent):
    Object3D(center, parent), origin(true), idxMaisProximo(-1),
    showInputLine(false), showTriQuad(true), meshTranslation(false), showScalarField(false), m_glTextureName(texName)
{
    setInputType(GL_TRIANGLES);
    buildObject();
}

TriQuadMesh::TriQuadMesh(const TriQuadMesh& tt): Object3D(tt), origin(false), idxMaisProximo(-1)
{
    m_glTextureName = tt.textureName();
    buildObject();
}

int TriQuadMesh::textureName()const
{
    return m_glTextureName;
}

bool TriQuadMesh::isEmpty()
{
    return che.sizeOfVertices() == 0;
}

void TriQuadMesh::clear()
{
    che.clear();
}

void TriQuadMesh::viewGrad(bool v)
{
    if(v)
        activeProgram = 1;
    else
        activeProgram = 0;
}

bool TriQuadMesh::buildMesh(CHEBuilder *builder)
{
    builder->build();
    che = builder->che();
    che.maxMimCalc();
    return !che.isEmpty();

//    che.clear();

//    Vertex v1(0,-0.3);
//    v1.quadric() = Quadric2D(4,0,0,0,6,0);

//    Vertex v2( 1, -0.3);
//    v2.quadric() = Quadric2D(0,0,0,0,1,0);

//    Vertex v3(0 , 1);
//    v3.quadric() = Quadric2D(0,0,0,0,1,0);//Quadric2D(1,0,0,1,0,1);

//    QVector<Vertex> v;
//    v.append(v1);v.append(v2);v.append(v3);
//    che.addVertices(v);
//    che.addTriangle(0,1,2);
//    return true;
}

Object3D* TriQuadMesh::copy() const
{
    return new TriQuadMesh(*this);
}

TriQuadMesh::~TriQuadMesh()
{
}

void TriQuadMesh::beforeTransformations(void)
{
    //if(origin)
    //    drawOrigin();
}

void TriQuadMesh::clearDrawPoints()
{
    dp.clear();
    dpL.clear();
    dpU.clear();
}

void TriQuadMesh::drawPoints(const QVector<QVector2D>& ps)
{
    if(ps.size() > 0)
        dp = ps;

    glColor4f(1.0, 0.0, 0.0, 1.0);
    glPointSize(4.0);
    glBegin(GL_POINTS);
    for(int i = 0; i < dp.size(); ++i)
        glVertex3f(dp[i].x(), dp[i].y(), 2.0);
    glEnd();
    glPointSize(1.0);
}

void TriQuadMesh::drawPoints1(const QVector<QVector2D>& ps)
{
    if(ps.size() > 0)
        dpU = ps;

    glColor4f(1.0, 1.0, .0, 0.5);
    glPointSize(5.0);
    glBegin(GL_POINTS);
    for(int i = 0; i < dpU.size(); ++i)
        glVertex3f(dpU[i].x(), dpU[i].y(), 2.0);
    glEnd();
    glPointSize(1.0);
}

void TriQuadMesh::drawPoints2(const QVector<QVector2D>& ps)
{
    if(ps.size() > 0)
        dpL = ps;

    glColor4f(1.0, 0.0, 1.0, 0.5);
    glPointSize(5.0);
    glBegin(GL_POINTS);
    for(int i = 0; i < dpL.size(); ++i)
        glVertex3f(dpL[i].x(), dpL[i].y(), 2.0);
    glEnd();
    glPointSize(1.0);
}

void TriQuadMesh::drawGeometry(void)
{
    if(showSketch)
    {
        drawPoints();
        drawPoints1();
        drawPoints2();
    }

    if(showInputLine)
    {
        glColor4f(1.0, 0.0, 0.0, 1.0);
        glPointSize(1.0);
        glBegin(GL_LINE_STRIP);
        for(int i = 0; i < dp.size(); ++i)
            glVertex3f(dp[i].x(), dp[i].y(), 2.0);
        glEnd();
    }


    if(showMesh)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glLineWidth(2);
        glColor3f(0.7, 0.7, 0.7);
        glBegin(GL_TRIANGLES);
        {
            for(int i = 0; i < che.sizeOfTriangles(); ++i)
            {
                for(int j = 0; j < 3; ++j)
                {
                    int k = che.vertexId(i,j);
                    glVertex3f(che.vertex( k ).x(), che.vertex( k ).y(), 2.0);
                }
            }
        }glEnd();
        glLineWidth(1);
    }

    if(showTriQuad)
    {

        SLGl3W::activeTexture();
        //glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, m_glTextureName);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        program[activeProgram].bind();
        program[activeProgram].setUniformValue(locationScalar[activeProgram],showScalarField);
        program[activeProgram].setUniformValue(locationMaior[activeProgram], che.getMaior());
        program[activeProgram].setUniformValue(locationMenor[activeProgram], che.getMenor());
        program[activeProgram].setUniformValue(locationTexture[activeProgram], 0);
        glBegin(GL_TRIANGLES);
        {
            for(int i = 0; i < che.sizeOfTriangles(); ++i)
            {
                int v0 = che.vertexId(i,0);
                int v1 = che.vertexId(i,1);
                int v2 = che.vertexId(i,2);
                QMatrix4x4 W = buildInv(i);
                Quadric2D Qx = che.vertex(v0).quadric()* W(0,0) + che.vertex(v1).quadric()* W(1,0) + che.vertex(v2).quadric()* W(2,0);
                Quadric2D Qy = che.vertex(v0).quadric()* W(0,1) + che.vertex(v1).quadric()* W(1,1) + che.vertex(v2).quadric()* W(2,1);

                program[activeProgram].setUniformValue(locationQxABC[activeProgram], Qx.abc());
                program[activeProgram].setUniformValue(locationQxDEF[activeProgram], Qx.def());
                program[activeProgram].setUniformValue(locationQyABC[activeProgram], Qy.abc());
                program[activeProgram].setUniformValue(locationQyDEF[activeProgram], Qy.def());
                for(int j = 0; j < 3; ++j)
                {
                    int k = che.vertexId(i,j);
                    program[activeProgram].setAttributeValue(locationABC[activeProgram], che.vertex(k).quadric().abc());
                    program[activeProgram].setAttributeValue(locationDEF[activeProgram], che.vertex(k).quadric().def());

                    glVertex2f(che.vertex(k).x(),che.vertex(k).y());

                }
            }

        }glEnd();

        glBindTexture(GL_TEXTURE_2D,0);

        program[activeProgram].release();
    }

}

void TriQuadMesh::buildObject()
{
    showMesh = true;
    showSketch = true;
    showScalarField = true;

    QGLShader *vert = new QGLShader(QGLShader::Vertex  );
    QGLShader *vertG = new QGLShader(QGLShader::Vertex  );
    QGLShader *frag = new QGLShader(QGLShader::Fragment);
    QGLShader *fragG = new QGLShader(QGLShader::Fragment);
    //QGLShader *geom = new QGLShader(QGLShader::Geometry);

    vert->compileSourceFile(":/triquadVert");
    vert->log();

    vertG->compileSourceFile(":/triquadVert");
    vertG->log();

    frag->compileSourceFile(":/triquadFrag");
    frag->log();

    fragG->compileSourceFile(":/gradVisFrag");
    fragG->log();

    //geom->compileSourceFile(":/gradVisGeom");
    //geom->log();

    program[0].addShader(vert);
    program[0].addShader(frag);

    program[1].addShader(vertG);
    program[1].addShader(fragG);
    //program[1].addShader(geom);

    //program[1].setGeometryInputType(GL_TRIANGLES);
    //program[1].setGeometryOutputType(GL_TRIANGLES);
    //program[1].setGeometryOutputVertexCount(3);

    program[0].link();
    program[0].log();

    program[1].link();
    program[1].log();

    locationABC[0] = program[0].attributeLocation("abc");
    locationDEF[0] = program[0].attributeLocation("def");
    locationScalar[0] = program[0].uniformLocation("showScalar");
    locationMenor[0] = program[0].uniformLocation("menor");
    locationMaior[0] = program[0].uniformLocation("maior");
    locationTexture[0] = program[0].uniformLocation("sampler2d0");

    locationABC[1] = program[1].attributeLocation("abc");
    locationDEF[1] = program[1].attributeLocation("def");
    locationScalar[1] = program[1].uniformLocation("showScalar");
    locationMenor[1] = program[1].uniformLocation("menor");
    locationMaior[1] = program[1].uniformLocation("maior");
    locationTexture[1] = program[1].uniformLocation("sampler2d0");
    locationQxABC[1] = program[1].uniformLocation("Qxabc");
    locationQxDEF[1] = program[1].uniformLocation("Qxdef");
    locationQyABC[1] = program[1].uniformLocation("Qyabc");
    locationQyDEF[1] = program[1].uniformLocation("Qydef");

    activeProgram = 0;

    SLGl3W::init();
}

bool TriQuadMesh::isProgramLinked()
{
    return program[activeProgram].isLinked();
}

void TriQuadMesh::viewMesh(bool v)
{
    showMesh = v;
}

void TriQuadMesh::viewSketch(bool v)
{
    showSketch = v;
}

void TriQuadMesh::drawOrigin()
{
    GLboolean isLighting;
    GLfloat color[4];
    glGetBooleanv(GL_LIGHTING,&isLighting);
    glGetFloatv(GL_CURRENT_COLOR,color);
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);

    glBegin(GL_LINES);
    glColor3f(1.0, 0, 0);
    glVertex3i(0, 0, 0);
    glVertex3i(1, 0, 0);

    glColor3f(0, 1.0, 0);
    glVertex3i(0, 0, 0);
    glVertex3i(0, 1, 0);
    glEnd();

    glEnable(GL_TEXTURE_2D);
    glColor4fv(color);
    if(isLighting)
        glEnable(GL_LIGHTING);
}

void TriQuadMesh::changeOrigin(bool v)
{
    origin = v;
}

QMatrix4x4 TriQuadMesh::glGetMatrix(GLenum fetchType)
{
    QMatrix4x4 ret;
    GLfloat mat[16];
    glGetFloatv(fetchType, mat);
    qreal *m = ret.data();
    for (int index = 0; index < 16; ++index)
        m[index] = mat[index];

    return ret;
}
void TriQuadMesh::afterTransformations(void)
{
    glGetIntegerv(GL_VIEWPORT, vwp);

    QMatrix4x4 mvm = glGetMatrix(GL_MODELVIEW_MATRIX);
    QMatrix4x4 pjm = glGetMatrix(GL_PROJECTION_MATRIX);

    mvpi = (pjm*mvm).inverted();
}

QVector4D TriQuadMesh::unproject(const QPoint& p)
{
    QVector4D v(p);
    v.setW( 1.0 );
    v.setX(2.0*(v.x() - vwp[0]));
    v.setY(2.0*(v.y() - vwp[1]));
    v.setX(v.x()/vwp[2] -1.0);
    v.setY(-(v.y()/vwp[3] -1.0));

    v = mvpi*v;
    v /= v.w();
    v.setZ(1.0);
    return v;
}

QVector4D TriQuadMesh::unproject(const QPointF& p)
{
    QVector4D v(p);
    v.setW( 1.0 );
    v.setX(2.0*(v.x() - vwp[0]));
    v.setY(2.0*(v.y() - vwp[1]));
    v.setX(v.x()/vwp[2] -1.0);
    v.setY(-(v.y()/vwp[3] -1.0));

    v = mvpi*v;
    v /= v.w();
    v.setZ(1.0);
    return v;
}

QVector<QVector4D> TriQuadMesh::unproject(const QVector<QPoint>& v)
{
    QVector<QVector4D> ret;

    for(int i = 0; i < v.size(); ++i)
    {
        ret.append(unproject(v[i]));
    }

    return ret;
}

QVector<QVector4D> TriQuadMesh::unproject(const QVector<QPointF>& v)
{
    QVector<QVector4D> ret;

    for(int i = 0; i < v.size(); ++i)
    {
        ret.append(unproject(v[i]));
    }

    return ret;
}

void TriQuadMesh::move(const QPoint& ini, const QPoint& curr)
{
    static QVector2D currA;
    if(!meshTranslation)
    {
        if(idxMaisProximo < 0)
        {
            QVector4D v = unproject(ini);
            idxMaisProximo = che.mostClosedVertex(v.toVector2D());
            if(idxMaisProximo >= 0)
                maisProximo =  che.vertex(idxMaisProximo).toVector2D();
        }else
        {
            che.vertex(idxMaisProximo) = maisProximo + (unproject(curr)-unproject(ini)).toVector2D();
        }
    }else
    {
        if(idxMaisProximo < 0)
        {
            currA = unproject(ini).toVector2D();
            idxMaisProximo = 0;
        }else
        {
            QVector2D v = unproject(curr).toVector2D();
            che.addToAllVertices(v - currA);
            currA = v;
        }
    }

    che.maxMimCalc();

}

void TriQuadMesh::finish()
{
    idxMaisProximo = -1;
    che.maxMimCalc();
}

void TriQuadMesh::cancel()
{
    if(!meshTranslation)
    {
        if(idxMaisProximo >= 0)
            che.vertex(idxMaisProximo) = maisProximo;
    }
    idxMaisProximo = -1;
    che.maxMimCalc();
}

int TriQuadMesh::configPoints( QVector<QVector4D>& pontos, QVector<QVector3D>& bary, QVector<int>& idx )
{
    QVector<int> remover;
    bary.clear();
    idx.clear();
    QVector<QMatrix4x4> invs(che.sizeOfTriangles());
    for(int i = 0; i < che.sizeOfTriangles(); ++i)
        invs[i] = buildInv(i);

    for(int j = 0; j < pontos.size(); ++j)
    {
        bool entrou = false;
        for(int i = 0; i < che.sizeOfTriangles(); ++i)
        {
            QVector4D v = invs[i] * pontos[j];
            if(v.x() >= 0.0 && v.x() <= 1.0 && v.y() >= 0.0 && v.y() <= 1.0 && v.z() >= 0.0 && v.z() <= 1.0)
            {
                bary.push_back(v.toVector3D());
                idx.push_back(i);
                entrou = true;
                break;
            }
        }
        if(!entrou)
            remover.append(j);
    }

    for(int i = remover.size()-1; i>=0 ; --i)
        pontos.remove(remover[i]);

    return pontos.size();
}

QVector3D TriQuadMesh::bary(const QVector4D& p, int *idx)
{
    for(int i = 0; i < che.sizeOfTriangles(); ++i)
    {
        QVector4D v = buildInv(i) * p;
        if(v.x() >= 0.0 && v.x() <= 1.0 && v.y() >= 0.0 && v.y() <= 1.0 && v.z() >= 0.0 && v.z() <= 1.0)
        {
            *idx = i;
            return v.toVector3D();
        }
    }
    *idx = -1.0;
    return QVector3D();
}

void TriQuadMesh::globalFitting_1layer(QVector<QVector4D>  pontos,  bool includeVertices)
{
    QVector<QVector3D> b;
    QVector<int> idx;
    QVector<QVector2D> pontos2D;
    QVector<QVector2D> pontos2DL;
    QVector<QVector2D> pontos2DU;
    clearDrawPoints();


    int np = configPoints(pontos, b, idx);
    int nq = che.sizeOfVertices();
    int nc = 5;
    if(np < nq*nc)
        return;

    Curve c(pontos.size());
    for(int i = 0; i < pontos.size(); ++i)
    {
        c.x(i) = pontos[i].x();
        c.y(i) = pontos[i].y();
    }
    c.compute_nosso_curv(1);

    float extraPoints = 0.0;
    if(includeVertices)
        extraPoints = che.sizeOfVertices();

    gsl_matrix * A = gsl_matrix_calloc (np + extraPoints, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np + extraPoints);

    for (int i = 0; i < np; ++i)
    {
        int tId = idx[i];
        double bar[3];
        bar[0] = b[i].x();
        bar[1] = b[i].y();
        bar[2] = b[i].z();

        QVector2D p = pontos[i].toVector2D();
        pontos2D.append(p);

        gsl_vector_set( B, i , 1);
        for (int j = 0; j < 3; ++j)
        {
            int vId = che.vertexId(tId, j);
            gsl_matrix_set (A, i , vId*nc + 0,     p.x()*p.x()*bar[j]); //x^2
            gsl_matrix_set (A, i , vId*nc + 1, 2.0*p.x()*p.y()*bar[j]); //2xy
            gsl_matrix_set (A, i , vId*nc + 2, 2.0*p.x()      *bar[j]); //2x
            gsl_matrix_set (A, i , vId*nc + 3,     p.y()*p.y()*bar[j]); //y^2
            gsl_matrix_set (A, i , vId*nc + 4, 2.0*p.y()      *bar[j]); //2y
        }
    }

    if(includeVertices)
    {
        QVector<float> distances(che.sizeOfVertices());
        QVector<int> pointsIndex(che.sizeOfVertices());
        calcDistances(pontos, &distances, &pointsIndex);
        for(int i = 0; i < che.sizeOfVertices(); ++i)
        {
            Vertex v = che.vertex(i);
            int pIdx = pointsIndex[i];
            float dist = distances[i] * SIGN(QVector2D::dotProduct(v.toVector2D()-pontos[pIdx].toVector2D(), QVector2D(c.nx(pIdx), c.ny(pIdx))));
            gsl_matrix_set (A,np + i, i*nc + 0,     v.x()*v.x()); //x^2
            gsl_matrix_set (A,np + i, i*nc + 1, 2.0*v.x()*v.y()); //2xy
            gsl_matrix_set (A,np + i, i*nc + 2, 2.0*v.x()      ); //2x
            gsl_matrix_set (A,np + i, i*nc + 3,     v.y()*v.y()); //y^2
            gsl_matrix_set (A,np + i, i*nc + 4, 2.0*v.y()      ); //2y
            gsl_vector_set (B,np + i, dist+1);
            if(dist > 0)
                pontos2DU.append(v.toVector2D());
            else
                pontos2DL.append(v.toVector2D());

        }
    }

    QVector<Quadric2D> qs = fittingGLOBAL(A,B, -1.0);

    gsl_matrix_free (A);
    gsl_vector_free (B);

    for(int i = 0; i < qs.size(); ++i)
    {
        che.vertex(i).quadric() = qs[i];
    }
    drawPoints(pontos2D);
    drawPoints1(pontos2DU);
    drawPoints2(pontos2DL);
}

void TriQuadMesh::globalFitting_2layers(QVector<QVector4D> pontos, float k)
{
    QVector<QVector3D> b;
    QVector<int> idx;
    QVector<QVector2D> pontos2DL;
    QVector<QVector2D> pontos2DU;
    clearDrawPoints();


    int np = configPoints(pontos, b, idx);
    int nq = che.sizeOfVertices();
    int nc = 5;
    if(np < nq*nc)
        return;

    Curve c(pontos.size());
    for(int i = 0; i < pontos.size(); ++i)
    {
        c.x(i) = pontos[i].x();
        c.y(i) = pontos[i].y();
    }
    c.compute_nosso_curv(1);

    gsl_matrix * A = gsl_matrix_calloc (np*2, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np*2);

    for (int i = 1; i < np-1; ++i)
    {
        double bar[3];
        QVector4D p[3];

        p[0] = pontos[i] + (QVector2D(c.nx(i), c.ny(i)).normalized()*k).toVector4D();
        pontos2DU.append(p[0].toVector2D());
        p[2] = pontos[i] - (QVector2D(c.nx(i), c.ny(i)).normalized()*k).toVector4D();
        pontos2DL.append(p[2].toVector2D());

        int tId;

        for(int k = 0; k < 3; ++k)
        {
            if(k == 1)
            {
                continue ;
                bar[0] = b[i].x();
                bar[1] = b[i].y();
                bar[2] = b[i].z();
                tId = idx[i];
            }
            else
            {
                int idx;
                QVector3D bp = bary(p[k], &idx);
                if(idx < 0)
                    continue;
                bar[0] = bp.x();
                bar[1] = bp.y();
                bar[2] = bp.z();
                tId = idx;
            }

            for (int j = 0; j < 3; ++j)
            {
                int vId = che.vertexId(tId, j);
                gsl_matrix_set (A, i*2 + k/2, vId*nc + 0,     p[k].x()*p[k].x()*bar[j]); //x^2
                gsl_matrix_set (A, i*2 + k/2, vId*nc + 1, 2.0*p[k].x()*p[k].y()*bar[j]); //2xy
                gsl_matrix_set (A, i*2 + k/2, vId*nc + 2, 2.0*p[k].x()         *bar[j]); //2x
                gsl_matrix_set (A, i*2 + k/2, vId*nc + 3,     p[k].y()*p[k].y()*bar[j]); //y^2
                gsl_matrix_set (A, i*2 + k/2, vId*nc + 4, 2.0*p[k].y()         *bar[j]); //2y
                gsl_vector_set( B, i*2 + k/2, 2-k);
            }
        }
    }

    QVector<Quadric2D> qs = fittingGLOBAL(A,B, -1.0);

    gsl_matrix_free (A);
    gsl_vector_free (B);

    for(int i = 0; i < qs.size(); ++i)
    {
        che.vertex(i).quadric() = qs[i];
    }
    //drawPoints(pontos2D);
    drawPoints1(pontos2DU);
    drawPoints2(pontos2DL);
}

void TriQuadMesh::globalFitting_3layers(QVector<QVector4D> pontos, float k, bool includeVertices)
{
    QVector<QVector3D> b;
    QVector<int> idx;
    QVector<QVector2D> pontos2D;
    QVector<QVector2D> pontos2DL;
    QVector<QVector2D> pontos2DU;
    clearDrawPoints();


    int np = configPoints(pontos, b, idx);
    int nq = che.sizeOfVertices();
    int nc = 5;
    if(np < nq*nc)
        return;

    Curve c(pontos.size());
    for(int i = 0; i < pontos.size(); ++i)
    {
        c.x(i) = pontos[i].x();
        c.y(i) = pontos[i].y();
    }
    c.compute_nosso_curv(1);

    float extraPoints = 0.0;
    if(includeVertices)
        extraPoints = che.sizeOfVertices();

    gsl_matrix * A = gsl_matrix_calloc (np*3 + extraPoints, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np*3 + extraPoints);

    for (int i = 1; i < np-1; ++i)
    {
        double bar[3];
        QVector4D p[3];

        p[0] = pontos[i] + (QVector2D(c.nx(i), c.ny(i)).normalized()*k).toVector4D();
        pontos2DU.append(p[0].toVector2D());
        p[1] = pontos[i];
        pontos2D.append(p[1].toVector2D());
        p[2] = pontos[i] - (QVector2D(c.nx(i), c.ny(i)).normalized()*k).toVector4D();
        pontos2DL.append(p[2].toVector2D());

        int tId;

        for(int k = 0; k < 3; ++k)
        {
            if(k == 1)
            {
                bar[0] = b[i].x();
                bar[1] = b[i].y();
                bar[2] = b[i].z();
                tId = idx[i];
            }
            else
            {
                int idx;
                QVector3D bp = bary(p[k], &idx);
                if(idx < 0)
                    continue;
                bar[0] = bp.x();
                bar[1] = bp.y();
                bar[2] = bp.z();
                tId = idx;
            }

            for (int j = 0; j < 3; ++j)
            {
                int vId = che.vertexId(tId, j);
                gsl_matrix_set (A, i*3 + k, vId*nc + 0,     p[k].x()*p[k].x()*bar[j]); //x^2
                gsl_matrix_set (A, i*3 + k, vId*nc + 1, 2.0*p[k].x()*p[k].y()*bar[j]); //2xy
                gsl_matrix_set (A, i*3 + k, vId*nc + 2, 2.0*p[k].x()         *bar[j]); //2x
                gsl_matrix_set (A, i*3 + k, vId*nc + 3,     p[k].y()*p[k].y()*bar[j]); //y^2
                gsl_matrix_set (A, i*3 + k, vId*nc + 4, 2.0*p[k].y()         *bar[j]); //2y
                gsl_vector_set( B, i*3 + k, 2-k);

            }
        }
    }

    QVector<float> distances(che.sizeOfVertices());
    QVector<int> pointsIndex(che.sizeOfVertices());
    calcDistances(pontos, &distances, &pointsIndex);

    if(includeVertices)
    {
        for(int i = 0; i < che.sizeOfVertices(); ++i)
        {
            Vertex v = che.vertex(i);
            int pIdx = pointsIndex[i];
            float dist = distances[i] * SIGN(QVector2D::dotProduct(v.toVector2D()-pontos[pIdx].toVector2D(), QVector2D(c.nx(pIdx), c.ny(pIdx))));
            gsl_matrix_set (A,np*3 + i, i*nc + 0,     v.x()*v.x()); //x^2
            gsl_matrix_set (A,np*3 + i, i*nc + 1, 2.0*v.x()*v.y()); //2xy
            gsl_matrix_set (A,np*3 + i, i*nc + 2, 2.0*v.x()      ); //2x
            gsl_matrix_set (A,np*3 + i, i*nc + 3,     v.y()*v.y()); //y^2
            gsl_matrix_set (A,np*3 + i, i*nc + 4, 2.0*v.y()      ); //2y
            gsl_vector_set (B,np*3 + i, dist+1);
            if(dist > 0)
                pontos2DU.append(v.toVector2D());
            else
                pontos2DL.append(v.toVector2D());

        }
    }

    QVector<Quadric2D> qs = fittingGLOBAL(A,B, -1.0);

    gsl_matrix_free (A);
    gsl_vector_free (B);

    for(int i = 0; i < qs.size(); ++i)
    {
        che.vertex(i).quadric() = qs[i];
    }
    drawPoints(pontos2D);
    drawPoints1(pontos2DU);
    drawPoints2(pontos2DL);
}

void TriQuadMesh::globalFitting_5layers(QVector<QVector4D> pontos, float k)
{
    QVector<QVector3D> b;
    QVector<int> idx;
    QVector<QVector2D> pontos2D;
    QVector<QVector2D> pontos2DL;
    QVector<QVector2D> pontos2DU;
    clearDrawPoints();

    int np = configPoints(pontos, b, idx);
    int nq = che.sizeOfVertices();
    int nc = 5;
    if(np < nq*nc)
        return;

    Curve c(pontos.size());
    for(int i = 0; i < pontos.size(); ++i)
    {
        c.x(i) = pontos[i].x();
        c.y(i) = pontos[i].y();
    }
    c.compute_nosso_curv(1);

    gsl_matrix * A = gsl_matrix_calloc (np*5, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np*5);

    for (int i = 1; i < np-1; ++i)
    {
        double bar[3];
        QVector4D p[5];

        p[0] = pontos[i] + (QVector2D(c.nx(i), c.ny(i)).normalized()*k*2).toVector4D();
        pontos2DU.append(p[0].toVector2D());
        p[1] = pontos[i] + (QVector2D(c.nx(i), c.ny(i)).normalized()*k).toVector4D();
        pontos2DU.append(p[1].toVector2D());
        p[2] = pontos[i];
        pontos2D.append(p[2].toVector2D());
        p[3] = pontos[i] - (QVector2D(c.nx(i), c.ny(i)).normalized()*k).toVector4D();
        pontos2DL.append(p[3].toVector2D());
        p[4] = pontos[i] - (QVector2D(c.nx(i), c.ny(i)).normalized()*k*2).toVector4D();
        pontos2DL.append(p[4].toVector2D());

        int tId;

        for(int k = 0; k < 5; ++k)
        {
            if(k == 2)
            {
                bar[0] = b[i].x();
                bar[1] = b[i].y();
                bar[2] = b[i].z();
                tId = idx[i];
            }
            else
            {
                int idx;
                QVector3D bp = bary(p[k], &idx);
                if(idx < 0)
                    continue;
                bar[0] = bp.x();
                bar[1] = bp.y();
                bar[2] = bp.z();
                tId = idx;
            }

            for (int j = 0; j < 3; ++j)
            {
                int vId = che.vertexId(tId,j);
                gsl_matrix_set (A, i*5 + k, vId*nc + 0,     p[k].x()*p[k].x()*bar[j]); //x^2
                gsl_matrix_set (A, i*5 + k, vId*nc + 1, 2.0*p[k].x()*p[k].y()*bar[j]); //2xy
                gsl_matrix_set (A, i*5 + k, vId*nc + 2, 2.0*p[k].x()         *bar[j]); //2x
                gsl_matrix_set (A, i*5 + k, vId*nc + 3,     p[k].y()*p[k].y()*bar[j]); //y^2
                gsl_matrix_set (A, i*5 + k, vId*nc + 4, 2.0*p[k].y()         *bar[j]); //2y
                gsl_vector_set( B, i*5 + k, 3-k);
            }
        }
    }

    QVector<Quadric2D> qs = fittingGLOBAL(A,B, -1.0);

    gsl_matrix_free (A);
    gsl_vector_free (B);

    for(int i = 0; i < qs.size(); ++i)
    {
        che.vertex(i).quadric() = qs[i];
    }
    drawPoints(pontos2D);
    drawPoints1(pontos2DU);
    drawPoints2(pontos2DL);
}

void TriQuadMesh::globalFitting_2layers_freef(QVector<QVector4D> pontos, float k)
{
    QVector<QVector3D> b;
    QVector<int> idx;
    QVector<QVector2D> pontos2DL;
    QVector<QVector2D> pontos2D;
    QVector<QVector2D> pontos2DU;
    clearDrawPoints();


    int np = configPoints(pontos, b, idx);
    int nq = che.sizeOfVertices();
    int nc = 6;
    if(np < nq*nc)
        return;

    Curve c(pontos.size());
    for(int i = 0; i < pontos.size(); ++i)
    {
        c.x(i) = pontos[i].x();
        c.y(i) = pontos[i].y();
    }
    c.compute_nosso_curv(1);

    gsl_matrix * A = gsl_matrix_calloc (np*2, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np*2);

    for (int i = 1; i < np-1; ++i)
    {
        double bar[3];
        QVector4D p[3];

        p[0] = pontos[i] + (QVector2D(c.nx(i), c.ny(i)).normalized()*k).toVector4D();
        pontos2DU.append(p[0].toVector2D());
        //pontos2D.append(pontos[i].toVector2D());
        p[2] = pontos[i] - (QVector2D(c.nx(i), c.ny(i)).normalized()*k).toVector4D();
        pontos2DL.append(p[2].toVector2D());

        int tId;

        for(int k = 0; k < 3; ++k)
        {
            if(k == 1)
            {
                continue ;
                bar[0] = b[i].x();
                bar[1] = b[i].y();
                bar[2] = b[i].z();
                tId = idx[i];
            }
            else
            {
                int idx;
                QVector3D bp = bary(p[k], &idx);
                if(idx < 0)
                    continue;
                bar[0] = bp.x();
                bar[1] = bp.y();
                bar[2] = bp.z();
                tId = idx;
            }

            for (int j = 0; j < 3; ++j)
            {
                int vId = che.vertexId(tId, j);
                gsl_matrix_set (A, i*2 + k/2, vId*nc + 0,     p[k].x()*p[k].x()*bar[j]); //x^2
                gsl_matrix_set (A, i*2 + k/2, vId*nc + 1, 2.0*p[k].x()*p[k].y()*bar[j]); //2xy
                gsl_matrix_set (A, i*2 + k/2, vId*nc + 2, 2.0*p[k].x()         *bar[j]); //2x
                gsl_matrix_set (A, i*2 + k/2, vId*nc + 3,     p[k].y()*p[k].y()*bar[j]); //y^2
                gsl_matrix_set (A, i*2 + k/2, vId*nc + 4, 2.0*p[k].y()         *bar[j]); //2y
                gsl_matrix_set (A, i*2 + k/2, vId*nc + 5,      1.0             *bar[j]); //1
                gsl_vector_set( B, i*2 + k/2, 1-k);
            }
        }
    }

    QVector<Quadric2D> qs = fittingGLOBAL_flivre(A,B);

    gsl_matrix_free (A);
    gsl_vector_free (B);

    for(int i = 0; i < qs.size(); ++i)
    {
        che.vertex(i).quadric() = qs[i];
    }
    drawPoints(pontos2D);
    drawPoints1(pontos2DU);
    drawPoints2(pontos2DL);
}

void TriQuadMesh::globalFittingG_3layers_freef(QVector<QVector4D> pontos, float kDistance, bool includeVertices)
{
    QVector<QVector3D> b;
    QVector<int> idx;
    QVector<QVector2D> pontos2DL;
    QVector<QVector2D> pontos2D;
    QVector<QVector2D> pontos2DU;
    clearDrawPoints();


    int np = configPoints(pontos, b, idx);
    int nq = che.sizeOfVertices();
    int nc = 6;
    if(np < nq*nc)
        return;

    Curve c(pontos.size());
    for(int i = 0; i < pontos.size(); ++i)
    {
        c.x(i) = pontos[i].x();
        c.y(i) = pontos[i].y();
    }
    c.compute_nosso_curv(1);

    float extraPoints = 0.0;
    if(includeVertices)
        extraPoints = che.sizeOfVertices();

    gsl_matrix * A = gsl_matrix_calloc (np*3 + extraPoints, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np*3 + extraPoints);

    for (int i = 1; i < np-1; ++i)
    {
        double bar[3];
        QVector4D p[3];

        p[0] = pontos[i] + (QVector2D(c.nx(i), c.ny(i)).normalized()*kDistance).toVector4D();
        p[1] = pontos[i];
        pontos2DU.append(p[0].toVector2D());
        pontos2D.append(pontos[i].toVector2D());
        p[2] = pontos[i] - (QVector2D(c.nx(i), c.ny(i)).normalized()*kDistance).toVector4D();
        pontos2DL.append(p[2].toVector2D());

        int tId;
        for(int k = 0; k < 3; ++k)
        {

            if(k == 1)
            {
                bar[0] = b[i].x();
                bar[1] = b[i].y();
                bar[2] = b[i].z();
                tId = idx[i];
            }
            else
            {
                int idx;
                QVector3D bp = bary(p[k], &idx);
                if(idx < 0)
                    continue;
                bar[0] = bp.x();
                bar[1] = bp.y();
                bar[2] = bp.z();
                tId = idx;
            }
            gsl_vector_set( B, i*3 + k, 1-k);
            for (int j = 0; j < 3; ++j)
            {
                int vId = che.vertexId(tId,j);
                gsl_matrix_set (A, i*3 + k, vId*nc + 0,     p[k].x()*p[k].x()*bar[j]); //x^2
                gsl_matrix_set (A, i*3 + k, vId*nc + 1, 2.0*p[k].x()*p[k].y()*bar[j]); //2xy
                gsl_matrix_set (A, i*3 + k, vId*nc + 2, 2.0*p[k].x()         *bar[j]); //2x
                gsl_matrix_set (A, i*3 + k, vId*nc + 3,     p[k].y()*p[k].y()*bar[j]); //y^2
                gsl_matrix_set (A, i*3 + k, vId*nc + 4, 2.0*p[k].y()         *bar[j]); //2y
                gsl_matrix_set (A, i*3 + k, vId*nc + 5,      1.0             *bar[j]); //1
            }
        }
    }

    QVector<float> distances(che.sizeOfVertices());
    QVector<int> pointsIndex(che.sizeOfVertices());
    calcDistances(pontos, &distances, &pointsIndex);

    if(includeVertices)
    {
        for(int i = 0; i < che.sizeOfVertices(); ++i)
        {
            Vertex v = che.vertex(i);
            int pIdx = pointsIndex[i];
            float dist = distances[i] * SIGN(QVector2D::dotProduct(v.toVector2D()-pontos[pIdx].toVector2D(), QVector2D(c.nx(pIdx), c.ny(pIdx))));
            gsl_matrix_set (A,np*3 + i, i*nc + 0,     v.x()*v.x()); //x^2
            gsl_matrix_set (A,np*3 + i, i*nc + 1, 2.0*v.x()*v.y()); //2xy
            gsl_matrix_set (A,np*3 + i, i*nc + 2, 2.0*v.x()      ); //2x
            gsl_matrix_set (A,np*3 + i, i*nc + 3,     v.y()*v.y()); //y^2
            gsl_matrix_set (A,np*3 + i, i*nc + 4, 2.0*v.y()      ); //2y
            gsl_matrix_set (A,np*3 + i, i*nc + 5, 1.0            ); //1
            gsl_vector_set (B,np*3 + i, dist+1);
            if(dist > 0)
                pontos2DU.append(v.toVector2D());
            else
                pontos2DL.append(v.toVector2D());

        }
    }

    QVector<Quadric2D> qs = fittingGLOBAL_flivre(A,B);

    gsl_matrix_free (A);
    gsl_vector_free (B);

    for(int i = 0; i < qs.size(); ++i)
    {
        che.vertex(i).quadric() = qs[i];
    }
    drawPoints(pontos2D);
    drawPoints1(pontos2DU);
    drawPoints2(pontos2DL);
}

void TriQuadMesh::globalFittingG_3layers_freef_withGrad(QVector<QVector4D> pontos, float kDistance, bool includeVertices)
{
    QVector<QVector3D> b;
    QVector<int> idx;
    QVector<QVector2D> pontos2DL;
    QVector<QVector2D> pontos2D;
    QVector<QVector2D> pontos2DU;
    clearDrawPoints();

    int np = configPoints(pontos, b, idx);
    QVector<QPair<Vertex, HalfEdge> > edges = che.getIntersections(pontos);
    int nq = che.sizeOfVertices();
    int nc = 6;
    int nr = edges.size();
    if(np+nr*2 < nq*nc)
        return;

    Curve c(pontos.size());
    for(int i = 0; i < pontos.size(); ++i)
    {
        c.x(i) = pontos[i].x();
        c.y(i) = pontos[i].y();
    }
    c.compute_nosso_curv(1);

    float extraPoints = 0.0;
    if(includeVertices)
        extraPoints = che.sizeOfVertices();

    gsl_matrix * A = gsl_matrix_calloc (np*3 + extraPoints + nr*2, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np*3 + extraPoints + nr*2);

    for (int i = 1; i < np-1; ++i)
    {
        double bar[3];
        QVector4D p[3];

        p[0] = pontos[i] + (QVector2D(c.nx(i), c.ny(i)).normalized()*kDistance).toVector4D();
        p[1] = pontos[i];
        pontos2DU.append(p[0].toVector2D());
        pontos2D.append(pontos[i].toVector2D());
        p[2] = pontos[i] - (QVector2D(c.nx(i), c.ny(i)).normalized()*kDistance).toVector4D();
        pontos2DL.append(p[2].toVector2D());

        int tId;
        for(int k = 0; k < 3; ++k)
        {

            if(k == 1)
            {
                bar[0] = b[i].x();
                bar[1] = b[i].y();
                bar[2] = b[i].z();
                tId = idx[i];
            }
            else
            {
                int idx;
                QVector3D bp = bary(p[k], &idx);
                if(idx < 0)
                    continue;
                bar[0] = bp.x();
                bar[1] = bp.y();
                bar[2] = bp.z();
                tId = idx;
            }
            gsl_vector_set( B, i*3 + k, 1-k);
            for (int j = 0; j < 3; ++j)
            {
                int vId = che.vertexId(tId,j);
                gsl_matrix_set (A, i*3 + k, vId*nc + 0,     p[k].x()*p[k].x()*bar[j]); //x^2
                gsl_matrix_set (A, i*3 + k, vId*nc + 1, 2.0*p[k].x()*p[k].y()*bar[j]); //2xy
                gsl_matrix_set (A, i*3 + k, vId*nc + 2, 2.0*p[k].x()         *bar[j]); //2x
                gsl_matrix_set (A, i*3 + k, vId*nc + 3,     p[k].y()*p[k].y()*bar[j]); //y^2
                gsl_matrix_set (A, i*3 + k, vId*nc + 4, 2.0*p[k].y()         *bar[j]); //2y
                gsl_matrix_set (A, i*3 + k, vId*nc + 5,      1.0             *bar[j]); //1
            }
        }
    }

    QVector<float> distances(che.sizeOfVertices());
    QVector<int> pointsIndex(che.sizeOfVertices());
    calcDistances(pontos, &distances, &pointsIndex);

    int offset = np*3;
    if(includeVertices)
    {
        for(int i = 0; i < che.sizeOfVertices(); ++i)
        {
            Vertex v = che.vertex(i);
            int pIdx = pointsIndex[i];
            float dist = distances[i] * SIGN(QVector2D::dotProduct(v.toVector2D()-pontos[pIdx].toVector2D(), QVector2D(c.nx(pIdx), c.ny(pIdx))));
            gsl_matrix_set (A,offset + i, i*nc + 0,     v.x()*v.x()); //x^2
            gsl_matrix_set (A,offset + i, i*nc + 1, 2.0*v.x()*v.y()); //2xy
            gsl_matrix_set (A,offset + i, i*nc + 2, 2.0*v.x()      ); //2x
            gsl_matrix_set (A,offset + i, i*nc + 3,     v.y()*v.y()); //y^2
            gsl_matrix_set (A,offset + i, i*nc + 4, 2.0*v.y()      ); //2y
            gsl_matrix_set (A,offset + i, i*nc + 5, 1.0            ); //1
            gsl_vector_set (B,offset + i, dist+1);
            if(dist > 0)
                pontos2DU.append(v.toVector2D());
            else
                pontos2DL.append(v.toVector2D());

        }
    }

    QVector<QMatrix4x4> invs(che.sizeOfTriangles());
    for(int i = 0; i < che.sizeOfTriangles(); ++i)
        invs[i] = buildInv(i);

    offset += extraPoints;
    pontos2D.clear();
    for(int i = 0; i < edges.size(); ++i)
    {
        int t1 = che.triangleIDFromHeId(edges[i].second.twinIndex());
        int t2 = che.triangleIDFromHeId(che.mesh()[edges[i].second.twinIndex()].twinIndex() );
        Vertex p = edges[i].first;

        QMatrix4x4 w1 = invs[t1];
        QMatrix4x4 w2 = invs[t2];

        for(int j = 0; j < 2; ++j) //equations | entry of grad vector
        {
                for(int k = 0; k < 3; ++k) // each vertex on triangle
                {
                    int vId1 = che.vertexId(t1, k);
                    int vId2 = che.vertexId(t2, k);

                    gsl_matrix_set (A, offset + i*2 + j, vId1*nc + 0,     p.x()*p.x()*w1(k,j)); //x^2
                    gsl_matrix_set (A, offset + i*2 + j, vId1*nc + 1, 2.0*p.x()*p.y()*w1(k,j)); //2xy
                    gsl_matrix_set (A, offset + i*2 + j, vId1*nc + 2, 2.0*p.x()      *w1(k,j)); //2x
                    gsl_matrix_set (A, offset + i*2 + j, vId1*nc + 3,     p.y()*p.y()*w1(k,j)); //y^2
                    gsl_matrix_set (A, offset + i*2 + j, vId1*nc + 4, 2.0*p.y()      *w1(k,j)); //2y
                    gsl_matrix_set (A, offset + i*2 + j, vId1*nc + 5,                 w1(k,j)); //1

                    gsl_matrix_set (A, offset + i*2 + j, vId2*nc + 0, -     p.x()*p.x()*w2(k,j)); //x^2
                    gsl_matrix_set (A, offset + i*2 + j, vId2*nc + 1, - 2.0*p.x()*p.y()*w2(k,j)); //2xy
                    gsl_matrix_set (A, offset + i*2 + j, vId2*nc + 2, - 2.0*p.x()      *w2(k,j)); //2x
                    gsl_matrix_set (A, offset + i*2 + j, vId2*nc + 3, -     p.y()*p.y()*w2(k,j)); //y^2
                    gsl_matrix_set (A, offset + i*2 + j, vId2*nc + 4, - 2.0*p.y()      *w2(k,j)); //2y
                    gsl_matrix_set (A, offset + i*2 + j, vId2*nc + 5, -                 w2(k,j)); //1

                }

                gsl_vector_set( B, offset + i*2 + j, 0.0);
        }

        pontos2D.append(p.toVector2D());
    }


    QVector<Quadric2D> qs = fittingGLOBAL_flivre(A,B);

    gsl_matrix_free (A);
    gsl_vector_free (B);

    for(int i = 0; i < qs.size(); ++i)
    {
        che.vertex(i).quadric() = qs[i];
    }
    drawPoints(pontos2D);
    drawPoints1(pontos2DU);
    drawPoints2(pontos2DL);
}

QVector<QVector2D> TriQuadMesh::pointsOnEdge(const HalfEdge& h, int c)const
{
    QVector<QVector2D> ret;
    int vId1 = h.vertexIndex();
    int vId2 = che.mesh()[h.twinIndex()].vertexIndex();

    Vertex v = che.vertex(vId1);
    QVector2D v1 = v.toVector2D();
    v = che.vertex(vId2);
    QVector2D v2 = v.toVector2D();

    float dt = 1.0/(c-1);
    for(int i = 0; i < c; ++i)
    {
        float t = i*dt;
        ret.append((1-t)*v1 + (t)*v2);
    }

    return ret;
}

void TriQuadMesh::globalFittingG_3layers_freef_kDistance(QVector<QVector4D> pontos, float kDistance, bool includeVertices)
{
    QVector<QVector3D> b;
    QVector<int> idx;
    QVector<QVector2D> pontos2DL;
    QVector<QVector2D> pontos2D;
    QVector<QVector2D> pontos2DU;
    clearDrawPoints();


    int np = configPoints(pontos, b, idx);
    int nq = che.sizeOfVertices();
    int nc = 6;
    if(np < nq*nc)
        return;

    Curve c(pontos.size());
    for(int i = 0; i < pontos.size(); ++i)
    {
        c.x(i) = pontos[i].x();
        c.y(i) = pontos[i].y();
    }
    c.compute_nosso_curv(1);

    float extraPoints = 0.0;
    if(includeVertices)
        extraPoints = che.sizeOfVertices();

    gsl_matrix * A = gsl_matrix_calloc (np*3 + extraPoints, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np*3 + extraPoints);

    for (int i = 1; i < np-1; ++i)
    {
        double bar[3];
        QVector4D p[3];

        p[0] = pontos[i] + (QVector2D(c.nx(i), c.ny(i)).normalized()*kDistance).toVector4D();
        p[1] = pontos[i];
        pontos2DU.append(p[0].toVector2D());
        pontos2D.append(pontos[i].toVector2D());
        p[2] = pontos[i] - (QVector2D(c.nx(i), c.ny(i)).normalized()*kDistance).toVector4D();
        pontos2DL.append(p[2].toVector2D());

        int tId;

        gsl_vector_set( B, i*3 + 0, kDistance);
        gsl_vector_set( B, i*3 + 1, 0.0);
        gsl_vector_set( B, i*3 + 2, -kDistance);
        for(int k = 0; k < 3; ++k)
        {

            if(k == 1)
            {
                bar[0] = b[i].x();
                bar[1] = b[i].y();
                bar[2] = b[i].z();
                tId = idx[i];
            }
            else
            {
                int idx;
                QVector3D bp = bary(p[k], &idx);
                if(idx < 0)
                    continue;
                bar[0] = bp.x();
                bar[1] = bp.y();
                bar[2] = bp.z();
                tId = idx;
            }


            for (int j = 0; j < 3; ++j)
            {
                int vId = che.vertexId(tId,j);
                gsl_matrix_set (A, i*3 + k, vId*nc + 0,     p[k].x()*p[k].x()*bar[j]); //x^2
                gsl_matrix_set (A, i*3 + k, vId*nc + 1, 2.0*p[k].x()*p[k].y()*bar[j]); //2xy
                gsl_matrix_set (A, i*3 + k, vId*nc + 2, 2.0*p[k].x()         *bar[j]); //2x
                gsl_matrix_set (A, i*3 + k, vId*nc + 3,     p[k].y()*p[k].y()*bar[j]); //y^2
                gsl_matrix_set (A, i*3 + k, vId*nc + 4, 2.0*p[k].y()         *bar[j]); //2y
                gsl_matrix_set (A, i*3 + k, vId*nc + 5,      1.0             *bar[j]); //1
            }
        }
    }

    QVector<float> distances(che.sizeOfVertices());
    QVector<int> pointsIndex(che.sizeOfVertices());
    calcDistances(pontos, &distances, &pointsIndex);

    if(includeVertices)
    {
        for(int i = 0; i < che.sizeOfVertices(); ++i)
        {
            Vertex v = che.vertex(i);
            int pIdx = pointsIndex[i];
            float dist = distances[i] * SIGN(QVector2D::dotProduct(v.toVector2D()-pontos[pIdx].toVector2D(), QVector2D(c.nx(pIdx), c.ny(pIdx))));
            gsl_matrix_set (A,np*3 + i, i*nc + 0,     v.x()*v.x()); //x^2
            gsl_matrix_set (A,np*3 + i, i*nc + 1, 2.0*v.x()*v.y()); //2xy
            gsl_matrix_set (A,np*3 + i, i*nc + 2, 2.0*v.x()      ); //2x
            gsl_matrix_set (A,np*3 + i, i*nc + 3,     v.y()*v.y()); //y^2
            gsl_matrix_set (A,np*3 + i, i*nc + 4, 2.0*v.y()      ); //2y
            gsl_matrix_set (A,np*3 + i, i*nc + 5, 1.0            ); //1
            gsl_vector_set (B,np*3 + i, dist);
            if(dist > 0)
                pontos2DU.append(v.toVector2D());
            else
                pontos2DL.append(v.toVector2D());

        }
    }

    QVector<Quadric2D> qs = fittingGLOBAL_flivre(A,B);

    gsl_matrix_free (A);
    gsl_vector_free (B);

    for(int i = 0; i < qs.size(); ++i)
    {
        che.vertex(i).quadric() = qs[i];
    }
    drawPoints(pontos2D);
    drawPoints1(pontos2DU);
    drawPoints2(pontos2DL);
}


void TriQuadMesh::globalFittingG_1layers_freef(QVector<QVector4D> pontos)
{
    QVector<QVector3D> b;
    QVector<int> idx;
    QVector<QVector2D> pontos2DL;
    QVector<QVector2D> pontos2D;
    QVector<QVector2D> pontos2DU;
    clearDrawPoints();


    int np = configPoints(pontos, b, idx);
    int nq = che.sizeOfVertices();
    int nc = 6;
    if(np < nq*nc)
        return;

    Curve c(pontos.size());
    for(int i = 0; i < pontos.size(); ++i)
    {
        c.x(i) = pontos[i].x();
        c.y(i) = pontos[i].y();
    }
    c.compute_nosso_curv(1);

    float extraPoints = che.sizeOfVertices();

    gsl_matrix * A = gsl_matrix_calloc (np + extraPoints, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np + extraPoints);

    for (int i = 0; i < np; ++i)
    {
        double bar[3];
        QVector4D p;

        p = pontos[i];
        pontos2D.append(p.toVector2D());

        int tId;

        bar[0] = b[i].x();
        bar[1] = b[i].y();
        bar[2] = b[i].z();
        tId = idx[i];


        for (int j = 0; j < 3; ++j)
        {
            int vId = che.vertexId(tId,j);
            gsl_matrix_set (A, i , vId*nc + 0,     p.x()*p.x()*bar[j]); //x^2
            gsl_matrix_set (A, i , vId*nc + 1, 2.0*p.x()*p.y()*bar[j]); //2xy
            gsl_matrix_set (A, i , vId*nc + 2, 2.0*p.x()      *bar[j]); //2x
            gsl_matrix_set (A, i , vId*nc + 3,     p.y()*p.y()*bar[j]); //y^2
            gsl_matrix_set (A, i , vId*nc + 4, 2.0*p.y()      *bar[j]); //2y
            gsl_matrix_set (A, i , vId*nc + 5,      1.0       *bar[j]); //1
        }
        gsl_vector_set( B, i , 0.0);
    }

    QVector<float> distances(che.sizeOfVertices());
    QVector<int> pointsIndex(che.sizeOfVertices());
    calcDistances(pontos, &distances, &pointsIndex);

    double weight = 0.1 ;
    for(int i = 0; i < che.sizeOfVertices(); ++i)
    {
        Vertex v = che.vertex(i);
        int pIdx = pointsIndex[i];
        float dist = distances[i] * SIGN(QVector2D::dotProduct(v.toVector2D()-pontos[pIdx].toVector2D(), QVector2D(c.nx(pIdx), c.ny(pIdx))));
        gsl_matrix_set (A,np + i, i*nc + 0, weight  * (    v.x()*v.x()) ); //x^2
        gsl_matrix_set (A,np + i, i*nc + 1, weight  * (2.0*v.x()*v.y()) ); //2xy
        gsl_matrix_set (A,np + i, i*nc + 2, weight  * (2.0*v.x()      ) ); //2x
        gsl_matrix_set (A,np + i, i*nc + 3, weight  * (    v.y()*v.y()) ); //y^2
        gsl_matrix_set (A,np + i, i*nc + 4, weight  * (2.0*v.y()      ) ); //2y
        gsl_matrix_set (A,np + i, i*nc + 5, weight  * (1.0            ) ); //1
        gsl_vector_set (B,np + i, weight  * dist);
        if(dist > 0)
            pontos2DU.append(v.toVector2D());
        else
            pontos2DL.append(v.toVector2D());

    }


    QVector<Quadric2D> qs = fittingGLOBAL_flivre(A,B);

    gsl_matrix_free (A);
    gsl_vector_free (B);

    for(int i = 0; i < qs.size(); ++i)
    {
        che.vertex(i).quadric() = qs[i];
    }
    drawPoints(pontos2D);
    drawPoints1(pontos2DU);
    drawPoints2(pontos2DL);
}


void TriQuadMesh::calcDistances(const QVector<QVector4D> &pontos, QVector<float>* distances, QVector<int>* idxPoint)
{
    QVector<bool> map(che.sizeOfVertices());
    for(int i = 0; i < map.size(); ++i)
        map[i] = false;

    for(int i = 0; i < che.sizeOfVertices(); ++i)
    {
        for(int j = 1; j < pontos.size()-1; ++j)
        {
            Vertex vertex = che.vertex(i);
            float d = vertex.distanceTo(pontos[j].toVector2D());
            if(!map[i] || d < distances->at(i))
            {
                map[i] = true;
                distances->operator [](i) = d;
                idxPoint-> operator [](i) = j;
            }
        }
    }
}

QMatrix4x4 TriQuadMesh::buildInv(int triangleId)
{
    QMatrix4x4 m;
    for(int i = 0; i < 3; ++i)
    {
        m(0,i) = che.vertex(che.vertexId(triangleId,i)).x();
        m(1,i) = che.vertex(che.vertexId(triangleId,i)).y();
        m(2,i) = 1.0;
        m(3,i) = m(i,3) = 0.0;
    }
    m(3,3) = 1.0;
    return m.inverted();
}

void TriQuadMesh::addVertex(const QVector4D& newVertex)
{
    Vertex v(newVertex.x(),newVertex.y());
    che.addVertex(v);
}

void TriQuadMesh::joinVerticesAt(const QVector4D& controlPoint)
{
    Vertex p = controlPoint.toVector2D();

    che.joinVerticesAt(p);
}

void TriQuadMesh::deleteTriangleWith(const Vertex& v)
{
    int tId = che.findTriangleWith(v);
    if(tId >= 0)
        che.deleteTriangle(tId);
}

QDebug operator<< (QDebug d, const Quadric2D &model)
{
    d.nospace() << "Quadric[ ";
    d.nospace() << "A="<< model.a() << ", B=" << model.b() << ", C=" << model.c();
    d.nospace() << ", D="<< model.d() << ", E=" << model.e() << ", F=" << model.f() << " ]";
    d.space() <<"";
    return d;
}

void TriQuadMesh::viewScalarField(bool v)
{
    showScalarField = v;
}

bool TriQuadMesh::loadMesh(const QString &filename)
{
    QFile f(filename);

    if(!f.open(QIODevice::ReadOnly))
        return false;

    QTextStream s(f.readAll());
    f.close();

    che.clear();
    QVector<Vertex> vertices;

    int qtd;
    s >> qtd;
    for(int i = 0; i < qtd; ++i)
    {
        float x,y;
        s >> x >> y;
        vertices.append(Vertex(x,y));
    }
    che.addVertices(vertices);
    s >> qtd;
    for(int i = 0; i < qtd; ++i)
    {
        int v1, v2, v3;
        s >> v1 >> v2 >> v3;
        che.addTriangle(v1,v2,v3);
    }

    return true;
}

bool TriQuadMesh::saveMesh(const QString &filename) const
{
    QFile f(filename);

    if(!f.open(QIODevice::WriteOnly))
        return false;

    QTextStream s(&f);

    s << che.sizeOfVertices() << "\n";
    for(int i = 0; i < che.sizeOfVertices(); ++i)
        s << che.vertex(i).x() << " " << che.vertex(i).y() << "\n";

    s << che.sizeOfTriangles() << "\n";
    for(int i = 0; i < che.sizeOfTriangles(); ++i)
        s << che.vertexId(i,0) << " " << che.vertexId(i,1) << " " << che.vertexId(i,2) << "\n";

    f.close();
    return true;
}

void TriQuadMesh::configureRenderInputLine()
{
    showInputLine = true;
    showMesh = false;
    showSketch = false;
    showScalarField = false;
    showTriQuad = false;
}



void TriQuadMesh::configureRenderTriQuad()
{
    showInputLine = false;
    showMesh = false;
    showSketch = false;
    showScalarField = false;
    showTriQuad = true;
}

void TriQuadMesh::reconfigure(bool mesh,bool sketch, bool field)
{
    showInputLine = false;
    showTriQuad = true;
    showMesh = mesh;
    showSketch = sketch;
    showScalarField = field;
}

void TriQuadMesh::setMeshTranslation(bool v)
{
    meshTranslation = v;
}
