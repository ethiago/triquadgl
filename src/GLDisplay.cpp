#include <QDebug>
#include <QToolTip>
#include "GLDisplay.h"
#include <QMatrix4x4>
#ifdef __APPLE__
    #include <OpenGL/glu.h>
#else
    #include <GL/glu.h>
#endif

#define ZOOMSTEP 0.05

GLDisplay::GLDisplay(QWidget *parent) : QGLWidget(parent),
    rigthPressedPoint(NULLPOINT),
    leftPressedPoint(NULLPOINT), zoom(1.0)
{
    setStyleSheet("border: 2px solid black;");

    setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);
}

GLDisplay::~GLDisplay()
{
}


void GLDisplay::initializeGL()
{
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_TEXTURE_2D);

    glShadeModel( GL_SMOOTH );
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );	// specify implementation-specific hints

    glClearDepth( 1.0 );					// specify the clear value for the depth buffer
    glDepthFunc( GL_LEQUAL );
    glDisable(GL_CULL_FACE);
    //glEnable(GL_CULL_FACE);
    glDisable(GL_LIGHTING);
    glCullFace(GL_BACK);
}

void GLDisplay::resizeGL(int width, int height)
{
    glViewport(0, 0, width, height);
}

float GLDisplay::xDist(float aspect)
{
    if(aspect < 1)
        return (0.5)*zoom;
    else
        return (0.5*aspect)*zoom;
}

float GLDisplay::yDist(float aspect)
{
    if(aspect > 1)
        return (0.5)*zoom;
    else
        return (0.5*(1/aspect))*zoom;
}

void GLDisplay::paintGL()
{
    float aspect = width()*1.0/height();

    float xdist = xDist(aspect);
    float ydist = yDist(aspect);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glOrtho(-xdist*3.0, xdist*3.0, -ydist*3.0, +ydist*3.0, 5.0, 50);
    //glFrustum(-xdist, xdist, -ydist, +ydist, 5.0, 50);

    gluLookAt(0,0, 15, 0,0,0, 0, 1,0);

    glMatrixMode(GL_MODELVIEW);

    glClearColor( 0.0, 0.0, 0.0, 1.0);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();

    emit drawModel();

    glFlush();
}

void GLDisplay::mousePressEvent ( QMouseEvent * event )
{
    // Garante que somente 1 botao estara' pressionado e
    //   cancela a operacao anterior caso o outro seja
    //   pressionado
    if(event->button() == Qt::RightButton)
    {
        if(leftPressedPoint != NULLPOINT)
        {
            leftPressedPoint = NULLPOINT;
            rigthPressedPoint = NULLPOINT;
            emit mouseCancel();
        }else
        {
            rigthPressedPoint = event->pos();
        }
        event->accept();
    }else if(event->button() == Qt::LeftButton)
    {
        if(rigthPressedPoint != NULLPOINT)
        {
            leftPressedPoint = NULLPOINT;
            rigthPressedPoint = NULLPOINT;
            emit mouseCancel();
        }else
        {
            leftPressedPoint = event->pos();
        }
        event->accept();
    }else
        event->ignore();
}

void GLDisplay::mouseReleaseEvent ( QMouseEvent * event )
{
    if(event->button() == Qt::RightButton)
    {
        if(rigthPressedPoint != NULLPOINT)
        {
            emit mouseRigthFinish(rigthPressedPoint, event->pos());
            rigthPressedPoint = NULLPOINT;
        }
        event->accept();
    }else if(event->button() == Qt::LeftButton)
    {
        if(leftPressedPoint != NULLPOINT)
        {
            emit mouseLefthFinish(leftPressedPoint, event->pos());
            leftPressedPoint = NULLPOINT;
        }
    }else
        event->ignore();
}

void GLDisplay::mouseMoveEvent(QMouseEvent *event)
{
    if(event->buttons() == Qt::RightButton)
    {
        if(rigthPressedPoint != NULLPOINT)
        {
            emit mouseRigthMove(rigthPressedPoint, event->pos());
        }
        event->accept();
    }else if(event->buttons() == Qt::LeftButton)
    {
        if(leftPressedPoint != NULLPOINT)
        {
            emit mouseLeftMove(leftPressedPoint, event->pos());
        }
        event->accept();
    }else
        event->ignore();
}

void GLDisplay::wheelEvent(QWheelEvent * event)
{
    int nsteps = event->delta() / (8*15);

    if (event->orientation() == Qt::Vertical) {
        zoom += ZOOMSTEP * nsteps;
        event->accept();
        updateGL();
    }else event->ignore();
 }

void GLDisplay::mouseDoubleClickEvent(QMouseEvent *event)
{
    if(event->buttons() == Qt::LeftButton)
    {
        emit mouseDoubleClickLeft(event->pos());
        event->accept();
        return;
    }

    event->ignore();
}
