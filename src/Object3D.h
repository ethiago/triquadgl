#ifndef OBJECT3D_H
#define OBJECT3D_H

#include <QObject>
#include <QVector3D>
#include <QQuaternion>
#include <QMatrix4x4>
#include <QVector>
#include <QtOpenGL>

class Object3D : public QObject
{
    Q_OBJECT

    QVector3D m_center;
    QVector3D m_translation;
    QMatrix4x4 m_rotations;
    QQuaternion m_interactiveQuartenion;

    QMatrix4x4 m_modelView;
    QMatrix4x4 m_projection;

    QSize m_windowSize;

    //for glu objects
    int m_slices;
    int m_stacks;

    int m_texture;

    int m_inputType;

    virtual void beforeTransformations(void) = 0;
    virtual void afterTransformations(void) = 0;
    virtual void drawGeometry(void) = 0;
    double degreeFromCos(double) const;

    QMatrix4x4 glGetMatrix(GLenum fetchType);

public:
    explicit Object3D(const QVector3D& center = QVector3D(0.0,0.0,0.0),
                      QObject *parent = 0);
    explicit Object3D(const Object3D& obj);


    int slices(void)const;
    int stacks(void)const;

    void setRotations(const QMatrix4x4&);
    const QMatrix4x4& rotations(void)const;
    void addRotation(const QQuaternion&);

    void setInteractiveQuartenion(const QQuaternion& quaternion);
    const QQuaternion& interactiveQuartenion(void)const;

    void setCenter(const QVector3D& center);
    const QVector3D& center(void)const;

    void setTraslation(const QVector3D& translation);
    const QVector3D& translation(void)const;

    void setTexture(int texture);
    int texture(void)const;

    int inputType(void)const;

    void draw(void);

    void cleanTransformations();

    QMatrix4x4 modelViewMatrix();
    QMatrix4x4 projectionMatrix();

    const QSize& windowSize();

protected:

    void setInputType(int input);
    void setStacks(int stacks);
    void setSlices(int slices);

private slots:
    void resizeWindow(const QSize& size);
};

#endif // OBJECT3D_H
