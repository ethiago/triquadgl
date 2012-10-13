#include "cubic.h"

Cubic::Cubic()
{
    bulild();
}

QVector2D& Cubic::operator [](int i)
{
    return m_v[i];
}

void Cubic::setConsts(float x3, float x2y, float xy2, float y3, float x2, float xy, float y2, float x, float y, float c )
{
    consts[0] = x3 ;
    consts[1] = x2y;
    consts[2] = xy2;
    consts[3] = y3 ;
    consts[4] = x2 ;
    consts[5] = xy ;
    consts[6] = y2 ;
    consts[7] = x  ;
    consts[8] = y  ;
    consts[9] = c  ;
}

void Cubic::bulild()
{
    m_v[0] = QVector2D( -1.0, -1.0);
    m_v[1] = QVector2D(  1.0, -1.0);
    m_v[2] = QVector2D(  1.0,  1.0);

    QGLShader *vert = new QGLShader(QGLShader::Vertex  );
    QGLShader *frag = new QGLShader(QGLShader::Fragment);
    qDebug() << vert->compileSourceFile(":/cubic.vert");
    qDebug() << vert->log();
    qDebug() << frag->compileSourceFile(":/cubic.frag");
    qDebug() << frag->log();

    program.addShader(vert);
    program.addShader(frag);

    qDebug() << program.link();
    qDebug() << program.log();

    constsLocation = program.uniformLocation("consts");
}

void Cubic::beforeTransformations(void)
{
}

void Cubic::afterTransformations(void)
{
}

void Cubic::drawGeometry(void)
{
    drawOrigin();
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3f(0.3, 0.3, 0.3);
    glBegin(GL_TRIANGLES);
    {
        for(int j = 0; j < 3; ++j)
        {
            glVertex2fv( reinterpret_cast<const GLfloat *>( &m_v[j] ) );
        }
    }glEnd();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    program.bind();
    program.setUniformValueArray(constsLocation, consts, 10, 1);
    glBegin(GL_TRIANGLES);
    {
        for(int j = 0; j < 3; ++j)
        {
            glVertex2fv(reinterpret_cast<const GLfloat *>(&m_v[j]));
        }
    }glEnd();

    program.release();
}

float * Cubic::getConsts()
{
    return consts;
}
