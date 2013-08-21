QT       += core gui opengl

TARGET = TriQuadGL
TEMPLATE = app

win32 {
    LIBS += -L"C:\Program Files\GnuWin32\lib" -lgsl -lgslcblas
    INCLUDEPATH += "C:\Program Files\GnuWin32\include"
}

unix{
    LIBS += -L/opt/local/lib/ -lgsl -lgslcblas -lGLU
    INCLUDEPATH += /opt/local/include
}

CONFIG-=app_bundle

HEADERS  += src/curveN/catmullrom.hpp \
        src/curveN/curven.hpp \
        src/curveN/vectorn.hpp \
    src/mainwindow.h\
    src/GLDisplay.h \
    src/rendercontroller.h \
    src/Object3D.h \
    src/sketchcontroller.h \
    src/triquadmesh.h \
    src/fitting.h \
    src/Curve.h \
    src/Curvature.h \
    src/quadric2d.h \
    src/vertex.h \
    src/halfedge.h \
    src/compacthalfedge.h \
    src/chebuilder.h \
    src/chebuilderregulargrid.h \
    src/chebuilderdefault.h \
    src/chebuilderquadtreefrompointcloud.h \
    src/chebuilderregulargridfrompointcloud.h \
    src/fastmarching.h

SOURCES  += src/main.cpp \
        src/curveN/catmullrom.hpp \
        src/curveN/curven.hpp \
        src/curveN/vectorn.hpp \
    src/mainwindow.cpp\
    src/GLDisplay.cpp \
    src/rendercontroller.cpp \
    src/Object3D.cpp \
    src/sketchcontroller.cpp \
    src/triquadmesh.cpp \
    src/fitting.cpp \
    src/Curve.cpp \
    src/quadric2d.cpp \
    src/vertex.cpp \
    src/halfedge.cpp \
    src/compacthalfedge.cpp \
    src/chebuilder.cpp \
    src/chebuilderregulargrid.cpp \
    src/chebuilderdefault.cpp \
    src/chebuilderquadtreefrompointcloud.cpp \
    src/chebuilderregulargridfrompointcloud.cpp \
    src/fastmarching.cpp


FORMS    += src/mainwindow.ui

RESOURCES += \
    Shaders.qrc

OTHER_FILES += shaders/TriQuad.vert \
        shaders/TriQuad.frag
