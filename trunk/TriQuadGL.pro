QT       += core gui opengl

TARGET = TriQuadGL
TEMPLATE = app

INCLUDEPATH += /opt/local/include

unix{
    LIBS += -L/opt/local/lib/ -lgsl -lgslcblas -lGLU
    #LIBS += -lGLU -lglut
}

SOURCES += src/main.cpp\
    src/mainwindow.cpp\
    src/GLDisplay.cpp \
    src/rendercontroller.cpp \
    src/Object3D.cpp \
    src/sketchcontroller.cpp \
    src/triquadmesh.cpp \
    src/fitting.cpp


HEADERS  += src/mainwindow.h\
    src/GLDisplay.h \
    src/rendercontroller.h \
    src/Object3D.h \
    src/sketchcontroller.h \
    src/triquadmesh.h \
    src/fitting.h

FORMS    += src/mainwindow.ui

RESOURCES += \
    Shaders.qrc

OTHER_FILES += shaders/TriQuad.vert \
        shaders/TriQuad.frag
