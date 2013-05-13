#include "triquadmesh.h"
#include <QGLShader>
#include <QVector3D>
#include <QtOpenGL>
#include "fitting.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "Curve.h"

TriQuadMesh::TriQuadMesh(const QVector3D& center, QObject *parent):
    Object3D(center, parent), origin(true), idxMaisProximo(-1)
{
    setInputType(GL_TRIANGLES);
    buildObject();
}

TriQuadMesh::TriQuadMesh(const TriQuadMesh& tt): Object3D(tt), origin(false), idxMaisProximo(-1)
{
    buildObject();
}

bool TriQuadMesh::isEmpty()
{
    return che.sizeOfVertices() == 0;
}

void TriQuadMesh::clear()
{
    che.clear();
}

bool TriQuadMesh::buildMesh(CHEBuilder *builder)
{
    builder->build();
    che = builder->che();
    return !che.isEmpty();
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

    glColor3f(1.0, 0.0, 0.0);
    glPointSize(5.0);
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

    glColor3f(1.0, 1.0, .0);
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

    glColor3f(1.0, 0.0, 1.0);
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

    if(showMesh)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glColor3f(0.3, 0.3, 0.3);
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
    }
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    program.bind();
    program.setUniformValue(locationScalar,showScalarField);
    glBegin(GL_TRIANGLES);
    {
        for(int i = 0; i < che.sizeOfTriangles(); ++i)
        {
            for(int j = 0; j < 3; ++j)
            {
                int k = che.vertexId(i,j);
                program.setAttributeValue(locationABC, che.vertex(k).quadric().abc());
                program.setAttributeValue(locationDEF, che.vertex(k).quadric().def());

                glVertex2f(che.vertex(k).x(),che.vertex(k).y());

            }
        }

    }glEnd();

    program.release();

}

void TriQuadMesh::buildObject()
{
    showMesh = true;
    showSketch = true;
    showScalarField = true;

    QGLShader *vert = new QGLShader(QGLShader::Vertex  );
    QGLShader *frag = new QGLShader(QGLShader::Fragment);
    qDebug() << vert->compileSourceFile(":/triquadVert");
    qDebug() << vert->log();
    qDebug() << frag->compileSourceFile(":/triquadFrag");
    qDebug() << frag->log();

    program.addShader(vert);
    program.addShader(frag);

    qDebug() << program.link();
    qDebug() << program.log();

    locationABC = program.attributeLocation("abc");
    locationDEF = program.attributeLocation("def");
    locationScalar = program.uniformLocation("showScalar");
}

bool TriQuadMesh::isProgramLinked()
{
    return program.isLinked();
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
    if(idxMaisProximo < 0)
    {
        QVector4D v = unproject(ini);
        idxMaisProximo = che.mostClosedVertex(v.toVector2D());
        if(idxMaisProximo >= 0)
            maisProximo =  che.vertex(idxMaisProximo).toVector2D();
    }

    if(idxMaisProximo >= 0)
        che.vertex(idxMaisProximo) = maisProximo + (unproject(curr)-unproject(ini)).toVector2D();
}

void TriQuadMesh::finish()
{
    idxMaisProximo = -1;
}

void TriQuadMesh::cancel()
{
    if(idxMaisProximo >= 0)
        che.vertex(idxMaisProximo) = maisProximo;
    idxMaisProximo = -1;
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


void TriQuadMesh::fittingG(QVector<QVector4D> & pontos)
{
    QVector<QVector3D> b;
    QVector<int> idx;
    QVector<QVector2D> pontos2D;

    int np = configPoints(pontos, b, idx);
    qDebug() << "PASSOU";
    int nq = che.sizeOfVertices();
    int nc = 5;
    if(np < nq*nc)
        return;

    gsl_matrix * A = gsl_matrix_calloc (np, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np);

    for (int i = 0; i < np; ++i)
    {
        double bary[3];
        QVector4D p = pontos[i];
        bary[0] = b[i].x();
        bary[1] = b[i].y();
        bary[2] = b[i].z();
        for (int j = 0; j < 3; ++j)
        {
            int vID = che.vertexId(idx[i],j);
            gsl_matrix_set (A, i, vID*nc    ,     p.x()*p.x()*bary[j]); //x^2
            gsl_matrix_set (A, i, vID*nc + 1, 2.0*p.x()*p.y()*bary[j]); //2xy
            gsl_matrix_set (A, i, vID*nc + 2, 2.0*p.x()      *bary[j]); //2x
            gsl_matrix_set (A, i, vID*nc + 3,     p.y()*p.y()*bary[j]); //y^2
            gsl_matrix_set (A, i, vID*nc + 4, 2.0*p.y()      *bary[j]); //2y
        }
        gsl_vector_set(B, i, 1.0);
        pontos2D.push_back(p.toVector2D());
    }

    QVector<Quadric2D> qs = fittingGLOBAL(A,B, -1.0);

    gsl_matrix_free (A);
    gsl_vector_free (B);

    for(int i = 0; i < qs.size(); ++i)
    {
        che.vertex(i).quadric() = qs[i];
    }
    drawPoints(pontos2D);
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

void TriQuadMesh::globalFitting_2layers(QVector<QVector4D> &pontos)
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
        float k = 0.2;//0.1/fabs(c.k(i));

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

void TriQuadMesh::globalFitting_2layers_freef(QVector<QVector4D> &pontos)
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
        float k = 0.2;//0.1/fabs(c.k(i));

        p[0] = pontos[i] + (QVector2D(c.nx(i), c.ny(i)).normalized()*k).toVector4D();
        pontos2DU.append(p[0].toVector2D());
        pontos2D.append(pontos[i].toVector2D());
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

void TriQuadMesh::globalFittingG_3layers_freef(QVector<QVector4D> &pontos)
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

    gsl_matrix * A = gsl_matrix_calloc (np*3, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np*3);

    for (int i = 1; i < np-1; ++i)
    {
        double bar[3];
        QVector4D p[3];
        float k = 0.2;//0.1/fabs(c.k(i));

        p[0] = pontos[i] + (QVector2D(c.nx(i), c.ny(i)).normalized()*k).toVector4D();
        p[1] = pontos[i];
        pontos2DU.append(p[0].toVector2D());
        pontos2D.append(pontos[i].toVector2D());
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
                int vId = che.vertexId(tId,j);
                gsl_matrix_set (A, i*3 + k, vId*nc + 0,     p[k].x()*p[k].x()*bar[j]); //x^2
                gsl_matrix_set (A, i*3 + k, vId*nc + 1, 2.0*p[k].x()*p[k].y()*bar[j]); //2xy
                gsl_matrix_set (A, i*3 + k, vId*nc + 2, 2.0*p[k].x()         *bar[j]); //2x
                gsl_matrix_set (A, i*3 + k, vId*nc + 3,     p[k].y()*p[k].y()*bar[j]); //y^2
                gsl_matrix_set (A, i*3 + k, vId*nc + 4, 2.0*p[k].y()         *bar[j]); //2y
                gsl_matrix_set (A, i*3 + k, vId*nc + 5,      1.0             *bar[j]); //1
                gsl_vector_set( B, i*3 + k, 1-k);
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

void TriQuadMesh::globalFitting_3layers(QVector<QVector4D> &pontos)
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

    gsl_matrix * A = gsl_matrix_calloc (np*3, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np*3);

    for (int i = 1; i < np-1; ++i)
    {
        double bar[3];
        QVector4D p[3];
        float k = 0.05;//0.1/fabs(c.k(i));

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

void TriQuadMesh::globalFitting_5layers(QVector<QVector4D> &pontos)
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
        float k = 0.05;//0.1/fabs(c.k(i));

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

void TriQuadMesh::fittingGG(QVector<QVector4D> & pontos)
{
    QVector<QVector2D> pontos2D;

    int np =  pontos.size();
    int nq = che.sizeOfVertices();
    int nt = che.sizeOfTriangles();

    int nc = 5;

    gsl_matrix * A = gsl_matrix_calloc (np*nt, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np*nt);

    QVector<QMatrix4x4> invs(che.sizeOfTriangles());
    for(int i = 0; i < che.sizeOfTriangles(); ++i)
        invs[i] = buildInv(i);

    for (int i = 0; i < np; ++i)
    {
        QVector4D p = pontos[i];
        pontos2D.append(p.toVector2D());
        double s = 0;
        for(int k = 0; k < nt; ++k)
        {
            QVector4D b = invs[k] * p;
            double bary[3];

            bary[0] = b.x();
            bary[1] = b.y();
            bary[2] = b.z();
            for (int j = 0; j < 3; ++j)
            {
                int vId = che.vertexId(k,j);
                gsl_matrix_set (A, nt*i + k, vId*nc    ,     p.x()*p.x()*bary[j]); //x^2
                gsl_matrix_set (A, nt*i + k, vId*nc + 1, 2.0*p.x()*p.y()*bary[j]); //2xy
                gsl_matrix_set (A, nt*i + k, vId*nc + 2, 2.0*p.x()      *bary[j]); //2x
                gsl_matrix_set (A, nt*i + k, vId*nc + 3,     p.y()*p.y()*bary[j]); //y^2
                gsl_matrix_set (A, nt*i + k, vId*nc + 4, 2.0*p.y()      *bary[j]); //2y
            }
            s+= b.x()+b.y()+b.z();
        }
        gsl_vector_set(B, i, s);
        pontos2D.push_back(p.toVector2D());
    }

    QVector<Quadric2D> qs = fittingGLOBAL(A,B, -1.0);

    gsl_matrix_free (A);
    gsl_vector_free (B);

    for(int i = 0; i < qs.size(); ++i)
    {
        che.vertex(i).quadric() = qs[i];
    }
    drawPoints(pontos2D);
}

void TriQuadMesh::globalFittingWithNormals(QVector<QVector4D> & pontos)
{
    QVector<QVector3D> b;
    QVector<int> idx;
    QVector<QVector2D> pontos2D;

    int np = configPoints(pontos, b, idx);
    np--;
    qDebug() << "PASSOU";
    int nq = che.sizeOfVertices();
    int nc = 5;
    int nLihasPerPoint = 3;
    if(np*nLihasPerPoint < nq*nc)
        return;

    gsl_matrix * A = gsl_matrix_calloc (np*nLihasPerPoint, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np*nLihasPerPoint);

    QVector<QMatrix4x4> invs(che.sizeOfTriangles());
    for(int i = 0; i < che.sizeOfTriangles(); ++i)
        invs[i] = buildInv(i);

    for (int i = 0; i < np; ++i)
    {
        double bary[3];
        QVector4D p = pontos[i];
        bary[0] = b[i].x();
        bary[1] = b[i].y();
        bary[2] = b[i].z();
        QMatrix4x4 W(invs[idx[i]]);
        QVector2D tang = pontos[i+1].toVector2D() - pontos[i].toVector2D();
        tang.normalize();
        for (int k = 0; k < 3; ++k)
        {
            int vId = che.vertexId(idx[i],k);
            //ponto
            gsl_matrix_set (A, i*nLihasPerPoint, vId*nc    ,     p.x()*p.x()*bary[k]); //x^2 * B_k
            gsl_matrix_set (A, i*nLihasPerPoint, vId*nc + 1, 2.0*p.x()*p.y()*bary[k]); //2xy * B_k
            gsl_matrix_set (A, i*nLihasPerPoint, vId*nc + 2, 2.0*p.x()      *bary[k]); //2x  * B_k
            gsl_matrix_set (A, i*nLihasPerPoint, vId*nc + 3,     p.y()*p.y()*bary[k]); //y^2 * B_k
            gsl_matrix_set (A, i*nLihasPerPoint, vId*nc + 4, 2.0*p.y()      *bary[k]); //2y  * B_k

            //dx
            gsl_matrix_set (A, i*nLihasPerPoint +1, vId*nc    ,     p.x()*p.x()*W(k,0) + 2.0*p.x()*      bary[k]); //x^2 * W(k,0) + 2x * bary[k]
            gsl_matrix_set (A, i*nLihasPerPoint +1, vId*nc + 1, 2.0*p.x()*p.y()*W(k,0) + 2.0*      p.y()*bary[k]); //2xy * W(k,0) + 2y * bary[k]
            gsl_matrix_set (A, i*nLihasPerPoint +1, vId*nc + 2, 2.0*p.x()      *W(k,0) + 2.0*            bary[k]); //2x  * W(k,0) + 2  * bary[k]
            gsl_matrix_set (A, i*nLihasPerPoint +1, vId*nc + 3,     p.y()*p.y()*W(k,0)                          ); //y^2 * W(k,0)
            gsl_matrix_set (A, i*nLihasPerPoint +1, vId*nc + 4, 2.0*p.y()      *W(k,0)                          ); //2y  * W(k,0)

            //dy
            gsl_matrix_set (A, i*nLihasPerPoint +2, vId*nc    ,     p.x()*p.x()*W(k,1)                          ); //x^2 * W(k,1)
            gsl_matrix_set (A, i*nLihasPerPoint +2, vId*nc + 1, 2.0*p.x()*p.y()*W(k,1) + 2.0*p.x()*      bary[k]); //2xy * W(k,1) + 2x * bary[k]
            gsl_matrix_set (A, i*nLihasPerPoint +2, vId*nc + 2, 2.0*p.x()      *W(k,1)                          ); //2x  * W(k,1)
            gsl_matrix_set (A, i*nLihasPerPoint +2, vId*nc + 3,     p.y()*p.y()*W(k,1) + 2.0*p.y()*      bary[k]); //y^2 * W(k,1) + 2y * bary[k]
            gsl_matrix_set (A, i*nLihasPerPoint +2, vId*nc + 4, 2.0*p.y()      *W(k,1) + 2.0*            bary[k]); //2y  * W(k,1) + 2  * bary[k]

        }
        gsl_vector_set(B, i*nLihasPerPoint, 1.0);
        gsl_vector_set(B, i*nLihasPerPoint +1, tang.x());
        gsl_vector_set(B, i*nLihasPerPoint +2, tang.y());
        pontos2D.push_back(p.toVector2D());
    }

    QVector<Quadric2D> qs = fittingGLOBAL(A,B, -1.0);

    gsl_matrix_free (A);
    gsl_vector_free (B);

    for(int i = 0; i < qs.size(); ++i)
    {
        che.vertex(i).quadric() = qs[i];
    }
    drawPoints(pontos2D);
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
