#include "triquadmesh.h"
#include <QGLShader>
#include <QVector3D>
#include <QtOpenGL>
#include "fitting.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
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

void TriQuadMesh::getTriangle(int i, QVector2D& v1, QVector2D& v2, QVector2D& v3)
{
    v1 = vertices[ triquads[i].idx[0] ];
    v2 = vertices[ triquads[i].idx[1] ];
    v3 = vertices[ triquads[i].idx[2] ];
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
    drawOrigin();
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

    quadrics.append(PARABOLA);
    quadrics.append(CIRCLE);
    quadrics.append(PARABOLA);

    NO no;
    no.idx[0] = 0; no.idx[1] = 1; no.idx[2] = 2;
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

void TriQuadMesh::fromCubic(float * consts)
{
    //consts ordem: x3 x2y xy2 y3 x2 xy y2 x y c
    QMatrix4x4 inv = buildInv(triquads[0]);
    qDebug() << inv;


}
