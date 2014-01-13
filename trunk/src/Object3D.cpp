#include "Object3D.h"

Object3D::Object3D(const QVector3D& center, QObject *parent)
    : QObject(parent), m_center(center), m_translation(QVector3D()),
      m_rotations(QMatrix4x4()) ,
      m_interactiveQuartenion(QQuaternion()),
      m_slices(50), m_stacks (50), m_texture(-1), m_inputType(GL_TRIANGLES)
{
}

Object3D::Object3D(const Object3D& obj)
    : QObject(obj.parent()), m_center(obj.center()),
      m_translation(obj.translation()), m_rotations(obj.rotations()),
      m_interactiveQuartenion(obj.interactiveQuartenion()),
      m_slices(obj.slices()), m_stacks(obj.stacks()),
      m_texture(obj.texture()), m_inputType(obj.inputType())

{
}

void Object3D::setRotations(const QMatrix4x4& rotations)
{
    m_rotations = rotations;
}

const QMatrix4x4& Object3D::rotations(void)const
{
    return m_rotations;
}

void Object3D::setInteractiveQuartenion(const QQuaternion& quaternion)
{
    m_interactiveQuartenion = quaternion;
}

const QQuaternion& Object3D::interactiveQuartenion(void)const
{
    return m_interactiveQuartenion;
}

void Object3D::setCenter(const QVector3D& center)
{
    m_center = center;
}

const QVector3D& Object3D::center(void)const
{
    return m_center;
}

void Object3D::setSlices(int slices)
{
    m_slices = slices;
}

int Object3D::slices(void)const
{
    return m_slices;
}

void Object3D::setStacks(int stacks)
{
    m_stacks = stacks;
}

int Object3D::stacks(void)const
{
    return m_stacks;
}
void Object3D::setTraslation(const QVector3D& translation)
{
    m_translation = translation;
}

const QVector3D& Object3D::translation(void)const
{
    return m_translation;
}

void Object3D::setTexture(int texture)
{
    m_texture = texture;
}

int Object3D::texture(void)const
{
    return m_texture;
}

void Object3D::addRotation(const QQuaternion& rotation)
{
    QMatrix4x4 m;
    m.rotate(degreeFromCos(rotation.scalar()), 0, 0, rotation.z());
    m_rotations = m*m_rotations;
}

void Object3D::draw(void)
{
    beforeTransformations();
    glMatrixMode(GL_MODELVIEW_MATRIX);
    glPushMatrix();

    glTranslatef(translation().x(), translation().y(), translation().z());
    glTranslatef(center().x(), center().y(), center().z());

    QQuaternion t = interactiveQuartenion();
    glRotatef(degreeFromCos(t.scalar()),
              0, 0, t.z());

    //glMultMatrixf(rotations().constData());
    glMultMatrixd(rotations().constData());

    m_modelView = glGetMatrix(GL_MODELVIEW_MATRIX);
    m_projection = glGetMatrix(GL_PROJECTION_MATRIX);

    afterTransformations();
    drawGeometry();

    glPopMatrix();
}

void Object3D::cleanTransformations()
{
    setTraslation(QVector3D());
    setCenter(QVector3D());

    setInteractiveQuartenion(QQuaternion());
    setRotations(QMatrix4x4());
}

int Object3D::inputType(void)const
{
    return m_inputType;
}

void Object3D::setInputType(int input)
{
    m_inputType = input;
}

double Object3D::degreeFromCos(double _cos)const
{
    return (acos(_cos)*180.0)/M_PI;
}

QMatrix4x4 Object3D::modelViewMatrix()
{
    return m_modelView;
}

QMatrix4x4 Object3D::projectionMatrix()
{
    return m_projection;
}

QMatrix4x4 Object3D::glGetMatrix(GLenum fetchType)
{
    QMatrix4x4 ret;
    GLfloat mat[16];
    glGetFloatv(fetchType, mat);
    //float *m = ret.data();
    qreal *m = ret.data();
    for (int index = 0; index < 16; ++index)
        m[index] = mat[index];

    return ret;
}

const QSize& Object3D::windowSize()
{
    return m_windowSize;
}

void Object3D::resizeWindow(const QSize& size)
{
    m_windowSize = size;
}
