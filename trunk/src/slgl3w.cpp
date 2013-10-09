#include "slgl3w.h"
#include "gl3w.h"

SLGl3W::SLGl3W()
{
}

void SLGl3W::init()
{
    gl3wInit();
}

void SLGl3W::activeTexture()
{
    glActiveTexture(GL_TEXTURE0);
}
