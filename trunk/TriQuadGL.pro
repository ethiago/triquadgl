QT       += core gui opengl

TARGET = TriQuadGL
TEMPLATE = app

unix{
    LIBS += -L/opt/local/lib/
    LIBS += -lGLU -lglut
}

SOURCES += src/main.cpp\
    src/mainwindow.cpp\
    src/GLDisplay.cpp \
    src/rendercontroller.cpp \
    src/Object3D.cpp \
    src/triquad.cpp


HEADERS  += src/mainwindow.h\
    src/GLDisplay.h \
    src/rendercontroller.h \
    src/Object3D.h \
    src/triquad.h

FORMS    += src/mainwindow.ui

RESOURCES += \
    Shaders.qrc

OTHER_FILES += shaders/TriQuad.vert \
        shaders/TriQuad.frag
