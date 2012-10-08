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

Object3D* TriQuadMesh::copy() const
{
    return new TriQuadMesh(*this);
}

TriQuadMesh::~TriQuadMesh()
{
}

void TriQuadMesh::beforeTransformations(void)
{
    if(origin)
        drawOrigin();
}

void TriQuadMesh::drawPoints(const QVector<QVector2D>& ps)
{
    static QVector<QVector2D> p;
    if(ps.size() > 0)
        p = ps;

    glColor3f(1.0,0.0,0.0);
    glBegin(GL_POINTS);
    for(int i = 0; i < p.size(); ++i)
        glVertex2fv(reinterpret_cast<const GLfloat *>(&p[i]));
    glEnd();
}

void TriQuadMesh::drawGeometry(void)
{
    drawPoints();

    if(showMesh)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glColor3f(0.3, 0.3, 0.3);
        glBegin(GL_TRIANGLES);
        {
            for(int i = 0; i < triquads.size(); ++i)
            {
                for(int j = 0; j < 3; ++j)
                {
                    int k = triquads[i].idx[j];
                    glVertex2fv(reinterpret_cast<const GLfloat *>(&(vertices[k])));
                }
            }
        }glEnd();
    }
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    program.bind();
    glBegin(GL_TRIANGLES);
    {
        for(int i = 0; i < triquads.size(); ++i)
        {
            for(int j = 0; j < 3; ++j)
            {
                int k = triquads[i].idx[j];
                program.setAttributeValue(locationABC, quadrics[k].a_b_c);
                program.setAttributeValue(locationDEF, quadrics[k].d_e_f);
                glVertex2fv(reinterpret_cast<const GLfloat *>(&(vertices[k])));

            }
        }

    }glEnd();

    program.release();

}

Quadric TriQuadMesh::makeQuadric(float x2, float y2, float xy, float x, float y, float c)
{
    Quadric q;
    q.a_b_c = QVector3D(x2,xy/2.0,x/2.0);
    q.d_e_f = QVector3D(y2, y/2.0,c);
    return q;
}

void TriQuadMesh::buildObject()
{
    showMesh = true;

    vertices.append(QVector2D( -1.0, -1.0));
    vertices.append(QVector2D(  1.0, -1.0));
    vertices.append(QVector2D(  1.0,  1.0));
    vertices.append(QVector2D( -1.0,  1.0));
    vertices.append(QVector2D(  2.0,  0.0));
    vertices.append(QVector2D( -2.0,  0.0));

    quadrics.append(PARABOLA);
    quadrics.append(CIRCLE);
    quadrics.append(PARABOLA);
    quadrics.append(CIRCLE2);
    quadrics.append(CIRCLE2);
    quadrics.append(PARABOLA);

    NO no;
    no.idx[0] = 0; no.idx[1] = 1; no.idx[2] = 2;
    triquads.append(no);
    no.idx[0] = 2; no.idx[1] = 3; no.idx[2] = 0;
    triquads.append(no);
    no.idx[0] = 4; no.idx[1] = 2; no.idx[2] = 1;
    triquads.append(no);
    no.idx[0] = 5; no.idx[1] = 0; no.idx[2] = 3;
    triquads.append(no);

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
}

void TriQuadMesh::setQuadric(int idx,  const Quadric& q)
{
    quadrics[idx] = q;
}

bool TriQuadMesh::isProgramLinked()
{
    return program.isLinked();
}

void TriQuadMesh::viewMesh(bool v)
{
    showMesh = v;
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

void TriQuadMesh::move(const QPoint& ini, const QPoint& curr)
{
    if(idxMaisProximo < 0)
    {
        idxMaisProximo = busca(unproject(ini));
        maisProximo = vertices[idxMaisProximo];
    }

    vertices[idxMaisProximo] = maisProximo + (unproject(curr)-unproject(ini)).toVector2D();
}

void TriQuadMesh::finish()
{
    idxMaisProximo = -1;
}

void TriQuadMesh::cancel()
{
    vertices[idxMaisProximo] = maisProximo;
    idxMaisProximo = -1;
}

int TriQuadMesh::busca(const QVector4D& p)
{
    int _idxMaisProximo = 0;
    float dist = (vertices[0]-p).length();

    for(int i = 1; i < vertices.size(); ++i)
    {
        float td = (vertices[i]-p).length();
        if(td < dist)
        {
            dist = td;
            _idxMaisProximo = i;
        }
    }
    return _idxMaisProximo;
}

int TriQuadMesh::configPoints(const QVector<QPoint> & in, QVector<QVector4D>& pontos, QVector<QVector3D>& bary, QVector<int>& idx )
{
    pontos.clear();
    bary.clear();
    idx.clear();
    QVector<QMatrix4x4> invs(triquads.size());
    for(int i = 0; i < triquads.size(); ++i)
        invs[i] = buildInv(triquads[i]);

    for(int j = 0; j < in.size(); ++j)
    {
        QVector4D p = unproject(in[j]);
        for(int i = 0; i < triquads.size(); ++i)
        {
            QVector4D v = invs[i] * p;
            if(v.x() >= 0.0 && v.x() <= 1.0 && v.y() >= 0.0 && v.y() <= 1.0 && v.z() >= 0.0 && v.z() <= 1.0)
            {
                pontos.push_back(p);
                bary.push_back(v.toVector3D());
                idx.push_back(i);
                break;
            }
        }
    }
    return pontos.size();
}

void TriQuadMesh::fittingG(const QVector<QPoint> & in)
{
    QVector<QVector4D> pontos;
    QVector<QVector3D> b;
    QVector<int> idx;
    QVector<QVector2D> pontos2D;

    int np = configPoints(in, pontos, b, idx);
    int nq = vertices.size();
    int nc = 5;
    if(np < nq*nc)
        return;

    gsl_matrix * A = gsl_matrix_calloc (np, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np);

    for (int i = 0; i < np; ++i)
    {
        double bary[3];
        QVector4D p = pontos[i];
        NO* no = &triquads[idx[i]];
        bary[0] = b[i].x();
        bary[1] = b[i].y();
        bary[2] = b[i].z();
        double s = bary[0] + bary[1] + bary[2];
        for (int j = 0; j < 3; ++j)
        {
            gsl_matrix_set (A, i, no->idx[j]*nc    ,     p.x()*p.x()*bary[j]); //x^2
            gsl_matrix_set (A, i, no->idx[j]*nc + 1, 2.0*p.x()*p.y()*bary[j]); //2xy
            gsl_matrix_set (A, i, no->idx[j]*nc + 2, 2.0*p.x()      *bary[j]); //2x
            gsl_matrix_set (A, i, no->idx[j]*nc + 3,     p.y()*p.y()*bary[j]); //y^2
            gsl_matrix_set (A, i, no->idx[j]*nc + 4, 2.0*p.y()      *bary[j]); //2y
        }
        gsl_vector_set(B, i, s);
        pontos2D.push_back(p.toVector2D());
    }

    QVector<Quadric> qs = fittingGLOBAL(A,B);

    gsl_matrix_free (A);
    gsl_vector_free (B);

    for(int i = 0; i < qs.size(); ++i)
    {
        quadrics[i] = qs[i];
    }
    drawPoints(pontos2D);
}

#define MAKEVEC(i) QVector4D()

void TriQuadMesh::fittingG2(const QVector<QPoint> & in)
{
    QVector<QVector4D> pontos;
    QVector<QVector3D> b;
    QVector<int> idx;
    QVector<QVector2D> pontos2D;

    int np = configPoints(in, pontos, b, idx);
    int nq = vertices.size();
    int nc = 5;
    if(np < nq*nc)
        return;

    Curve c(pontos.size());
    for(int i = 0; i < pontos.size(); ++i)
    {
        c.x(i) = pontos[i].x();
        c.y(i) = pontos[i].y();
    }
    c.compute_nosso_curv(3);

    gsl_matrix * A = gsl_matrix_calloc (np*3, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np*3);

    for (int i = 0; i < np; ++i)
    {
        double bary[3];
        QVector4D ponto = pontos[i];
        NO* no = &triquads[idx[i]];
        bary[0] = b[i].x();
        bary[1] = b[i].y();
        bary[2] = b[i].z();
        double s = bary[0] + bary[1] + bary[2];
        for (int j = 0; j < 3; ++j)
        {
            QVector4D p = ponto + (QVector2D(c.nx(i), c.ny(i)).normalized()*c.k(i)).toVector4D();
            gsl_matrix_set (A, i*3 + 0, no->idx[j]*nc + 0,     p.x()*p.x()*bary[j]); //x^2
            gsl_matrix_set (A, i*3 + 0, no->idx[j]*nc + 1, 2.0*p.x()*p.y()*bary[j]); //2xy
            gsl_matrix_set (A, i*3 + 0, no->idx[j]*nc + 2, 2.0*p.x()      *bary[j]); //2x
            gsl_matrix_set (A, i*3 + 0, no->idx[j]*nc + 3,     p.y()*p.y()*bary[j]); //y^2
            gsl_matrix_set (A, i*3 + 0, no->idx[j]*nc + 4, 2.0*p.y()      *bary[j]); //2y
            gsl_vector_set( B, i*3 + 0, 3);

            p = ponto ;
            gsl_matrix_set (A, i*3 + 1, no->idx[j]*nc + 0,     p.x()*p.x()*bary[j]); //x^2
            gsl_matrix_set (A, i*3 + 1, no->idx[j]*nc + 1, 2.0*p.x()*p.y()*bary[j]); //2xy
            gsl_matrix_set (A, i*3 + 1, no->idx[j]*nc + 2, 2.0*p.x()      *bary[j]); //2x
            gsl_matrix_set (A, i*3 + 1, no->idx[j]*nc + 3,     p.y()*p.y()*bary[j]); //y^2
            gsl_matrix_set (A, i*3 + 1, no->idx[j]*nc + 4, 2.0*p.y()      *bary[j]); //2y
            gsl_vector_set( B, i*3 + 1, 2);

            p = ponto - (QVector2D(c.nx(i), c.ny(i)).normalized()*c.k(i)).toVector4D();
            gsl_matrix_set (A, i*3 + 2, no->idx[j]*nc + 0,     p.x()*p.x()*bary[j]); //x^2
            gsl_matrix_set (A, i*3 + 2, no->idx[j]*nc + 1, 2.0*p.x()*p.y()*bary[j]); //2xy
            gsl_matrix_set (A, i*3 + 2, no->idx[j]*nc + 2, 2.0*p.x()      *bary[j]); //2x
            gsl_matrix_set (A, i*3 + 2, no->idx[j]*nc + 3,     p.y()*p.y()*bary[j]); //y^2
            gsl_matrix_set (A, i*3 + 2, no->idx[j]*nc + 4, 2.0*p.y()      *bary[j]); //2y
            gsl_vector_set( B, i*3 + 2, 1);
        }
        pontos2D.push_back(ponto.toVector2D());
    }

    QVector<Quadric> qs = fittingGLOBAL(A,B);

    gsl_matrix_free (A);
    gsl_vector_free (B);

    for(int i = 0; i < qs.size(); ++i)
    {
        quadrics[i] = qs[i];
    }
    drawPoints(pontos2D);
}

QMatrix4x4 TriQuadMesh::buildInv(NO& no)
{
    QMatrix4x4 m;
    for(int i = 0; i < 3; ++i)
    {
        m(0,i) = vertices[no.idx[i]].x();
        m(1,i) = vertices[no.idx[i]].y();
        m(2,i) = 1.0;
        m(3,i) = m(i,3) = 0.0;
    }
    m(3,3) = 1.0;
    return m.inverted();
}

QDebug operator<< (QDebug d, const Quadric &model)
{
    d.nospace() << "Quadric[ ";
    d.nospace() << "A="<< model.a_b_c.x() << ", B=" << model.a_b_c.y() << ", C=" << model.a_b_c.z();
    d.nospace() << ", D="<< model.d_e_f.x() << ", E=" << model.d_e_f.y() << ", F=" << model.d_e_f.z() << " ]";
    d.space() <<"";
    return d;
}
