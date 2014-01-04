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
    glPointSize(5.0);
    glBegin(GL_POINTS);
    for(int i = 0; i < p.size(); ++i)
        glVertex3f(p[i].x(), p[i].y(), 0.0);
    glEnd();
    glPointSize(1.0);
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

    QVector<QVector4D> in2;
    QVector2D p0 = QVector2D(-2.5, 0.0);
    QVector2D p1 = QVector2D(-1.0, 5.0);
    QVector2D p2 = QVector2D( 1.0,-5.0);
    QVector2D p3 = QVector2D( 2.5, 0.0);

    for(int i = 0; i < 200; ++i)
    {
        float t = i/200.0;
        QVector2D p = pow((1-t),3)*p0 + pow((1-t),2)*t*p1 + (1-t)*t*t*p2 + t*t*t*p3;
        in2.append(QVector4D(p, 1.0, 1.0));
    }

    for(int j = 0; j < in2.size(); ++j)
    {
        QVector4D p = in2[j];//unproject(in[j]);
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
        gsl_vector_set(B, i, 1.0);
        pontos2D.push_back(p.toVector2D());
    }

    QVector<Quadric> qs = fittingGLOBAL(A,B, -1.0);

    gsl_matrix_free (A);
    gsl_vector_free (B);

    for(int i = 0; i < qs.size(); ++i)
    {
        quadrics[i] = qs[i];
    }
    drawPoints(pontos2D);
}

QVector3D TriQuadMesh::bary(const QVector4D& p, int &idx)
{
    for(int i = 0; i < triquads.size(); ++i)
    {
        QVector4D v = buildInv(triquads[i]) * p;
        if(v.x() >= 0.0 && v.x() <= 1.0 && v.y() >= 0.0 && v.y() <= 1.0 && v.z() >= 0.0 && v.z() <= 1.0)
        {
            idx = i;
            return v.toVector3D();
        }
    }
    idx = -1.0;
    return QVector3D();
}

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
    c.compute_nosso_curv(1);

    gsl_matrix * A = gsl_matrix_calloc (np*3, nq*nc);
    gsl_vector * B = gsl_vector_calloc(np*3);

    for (int i = 0; i < np; ++i)
    {
        double bar[3];
        QVector4D p[3];
        float k = 0.2;//0.1/fabs(c.k(i));
        p[0] = pontos[i] + (QVector2D(c.nx(i), c.ny(i)).normalized()*k).toVector4D();
        p[1] = pontos[i];
        p[2] = pontos[i] - (QVector2D(c.nx(i), c.ny(i)).normalized()*k).toVector4D();
        NO* no;

        double s = bar[0] + bar[1] + bar[2];
        for(int k = 0; k < 3; ++k)
        {
            if(k == 1)
            {
                bar[0] = b[i].x();
                bar[1] = b[i].y();
                bar[2] = b[i].z();
                no = &triquads[idx[i]];
            }
            else
            {
                int idx;
                QVector3D bp = bary(p[k], idx);
                if(idx < 0)
                    continue;
                bar[0] = bp.x();
                bar[1] = bp.y();
                bar[2] = bp.z();
                no = &triquads[idx];
            }
            pontos2D.push_back(p[k].toVector2D());
            for (int j = 0; j < 3; ++j)
            {

                gsl_matrix_set (A, i*3 + k, no->idx[j]*nc + 0,     p[k].x()*p[k].x()*bar[j]); //x^2
                gsl_matrix_set (A, i*3 + k, no->idx[j]*nc + 1, 2.0*p[k].x()*p[k].y()*bar[j]); //2xy
                gsl_matrix_set (A, i*3 + k, no->idx[j]*nc + 2, 2.0*p[k].x()         *bar[j]); //2x
                gsl_matrix_set (A, i*3 + k, no->idx[j]*nc + 3,     p[k].y()*p[k].y()*bar[j]); //y^2
                gsl_matrix_set (A, i*3 + k, no->idx[j]*nc + 4, 2.0*p[k].y()         *bar[j]); //2y
                gsl_vector_set( B, i*3 + k, 2-k);
            }
        }
    }

    QVector<Quadric> qs = fittingGLOBAL(A,B, -1.0);

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


void TriQuadMesh::saveSketch(const QVector<QPoint>& points)
{
    static int cnt = 0;

    QFile f("scketch"+QString::number(cnt)+".txt");

    f.open(QIODevice::WriteOnly | QIODevice::Text);

    f.write((QString::number(points.size()) + "\n").toAscii());

    for(int i = 0;i < points.size(); ++i)
    {
        QVector4D p = unproject(points[i]);
        f.write((QString::number(p.x()) + " " + QString::number(p.y()) + "\n").toAscii());
    }
    f.close();
    ++cnt;
}

void TriQuadMesh::saveMesh()
{
    static int cnt = 0;

    QFile f("mesh"+QString::number(cnt)+".txt");

    f.open(QIODevice::WriteOnly | QIODevice::Text);

    f.write((QString::number(vertices.size()) + "\n").toAscii());
    for(int i = 0;i < vertices.size(); ++i)
    {
        QVector2D p = vertices[i];
        f.write((QString::number(p.x()) + " " + QString::number(p.y()) + "\n").toAscii());
    }
    f.write((QString::number(triquads.size()) + "\n").toAscii());
    for(int i = 0;i < triquads.size(); ++i)
    {
        int i0 = triquads[i].idx[0];
        int i1 = triquads[i].idx[1];
        int i2 = triquads[i].idx[2];
        f.write((QString::number(i0) + " " + QString::number(i1) + " " + QString::number(i2) + "\n").toAscii());
    }
    f.close();
    ++cnt;
}