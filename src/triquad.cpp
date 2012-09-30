#include "triquad.h"
#include <QGLShader>
#include <QVector3D>

#ifdef __APPLE__
    #include <OpenGL/gl.h>
    #include <OpenGL/glext.h>
#else
    #include <GL/gl.h>
#endif

TriQuad::TriQuad(const QVector3D& center, QObject *parent):
    Object3D(center, parent), origin(true), idxMaisProximo(-1)
{
    setInputType(GL_TRIANGLES);
    buildObject();
}

TriQuad::TriQuad(const TriQuad& tt): Object3D(tt), origin(false), idxMaisProximo(-1)
{
    buildObject();
}

Object3D* TriQuad::copy() const
{
    return new TriQuad(*this);
}

TriQuad::~TriQuad()
{
}

void TriQuad::beforeTransformations(void)
{
    if(origin)
        drawOrigin();
}

void TriQuad::drawGeometry(void)
{

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
                    glVertex4fv(reinterpret_cast<const GLfloat *>(&(vertices[k])));
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
                glVertex4fv(reinterpret_cast<const GLfloat *>(&(vertices[k])));

            }
        }

    }glEnd();

    program.release();

}

Quadric TriQuad::makeQuadric(float x2, float y2, float xy, float x, float y, float c)
{
    Quadric q;
    q.a_b_c = QVector3D(x2,xy,x);
    q.d_e_f = QVector3D(y2, y,c);
    return q;
}

void TriQuad::buildObject()
{
    showMesh = true;

    vertices.append(QVector4D( -1.0, -1.0, 0.0, 1.0));
    vertices.append(QVector4D(  1.0, -1.0, 0.0, 1.0));
    vertices.append(QVector4D(  1.0,  1.0, 0.0, 1.0));
    vertices.append(QVector4D( -1.0,  1.0, 0.0, 1.0));

    quadrics.append(PARABOLA);
    quadrics.append(CIRCLE);
    quadrics.append(PARABOLA);
    quadrics.append(CIRCLE2);

    NO no;
    no.idx[0] = 0; no.idx[1] = 1; no.idx[2] = 2;
    triquads.append(no);
    no.idx[0] = 2; no.idx[1] = 3; no.idx[2] = 0;
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

void TriQuad::setQuadric(int idx,  const Quadric& q)
{
    quadrics[idx] = q;
}

bool TriQuad::isProgramLinked()
{
    return program.isLinked();
}

void TriQuad::viewMesh(bool v)
{
    showMesh = v;
}

void TriQuad::drawOrigin()
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

void TriQuad::changeOrigin(bool v)
{
    origin = v;
}

QMatrix4x4 TriQuad::glGetMatrix(GLenum fetchType)
{
    QMatrix4x4 ret;
    GLfloat mat[16];
    glGetFloatv(fetchType, mat);
    qreal *m = ret.data();
    for (int index = 0; index < 16; ++index)
        m[index] = mat[index];

    return ret;
}
void TriQuad::afterTransformations(void)
{
    glGetIntegerv(GL_VIEWPORT, vwp);

    QMatrix4x4 mvm = glGetMatrix(GL_MODELVIEW_MATRIX);
    QMatrix4x4 pjm = glGetMatrix(GL_PROJECTION_MATRIX);

    mvpi = (pjm*mvm).inverted();
}

QVector4D TriQuad::unproject(const QPoint& p)
{
    QVector4D v(p);
    v.setW( 1.0 );
    v.setX(2.0*(v.x() - vwp[0]));
    v.setY(2.0*(v.y() - vwp[1]));
    v.setX(v.x()/vwp[2] -1.0);
    v.setY(-(v.y()/vwp[3] -1.0));

    return mvpi*v;
}

void TriQuad::move(const QPoint& ini, const QPoint& curr)
{
    if(idxMaisProximo < 0)
    {
        idxMaisProximo = busca(unproject(ini));
        maisProximo = vertices[idxMaisProximo];
    }

    vertices[idxMaisProximo] = maisProximo + (unproject(curr)-unproject(ini));
}

void TriQuad::finish()
{
    idxMaisProximo = -1;
}

void TriQuad::cancel()
{
    vertices[idxMaisProximo] = maisProximo;
    idxMaisProximo = -1;
}

int TriQuad::busca(const QVector4D& p)
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
