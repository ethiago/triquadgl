QT       += core gui opengl

TARGET = TriQuadGLFromCubic
TEMPLATE = app

INCLUDEPATH += /opt/local/include

win32 {
    LIBS += -L"C:\Program Files\GnuWin32\lib"
    INCLUDEPATH += "C:\Program Files\GnuWin32\include"
}

unix{
    LIBS += -L/opt/local/lib/
}

LIBS += -lgsl -lgslcblas #-lGLU

SOURCES += src/main.cpp\
    src/mainwindow.cpp\
    src/GLDisplay.cpp \
    src/rendercontroller.cpp \
    src/Object3D.cpp \
    src/sketchcontroller.cpp \
    src/triquadmesh.cpp \
    src/fitting.cpp \
    src/Curve.cpp \
    src/cubic.cpp


HEADERS  += src/mainwindow.h\
    src/GLDisplay.h \
    src/rendercontroller.h \
    src/Object3D.h \
    src/sketchcontroller.h \
    src/triquadmesh.h \
    src/fitting.h \
    src/Curve.h \
    src/Curvature.h \
    src/cubic.h

FORMS    += src/mainwindow.ui

RESOURCES += \
    Shaders.qrc

OTHER_FILES += shaders/TriQuad.vert \
        shaders/TriQuad.frag \
        shaders/cubic.vert \
        shaders/cubic.frag \
