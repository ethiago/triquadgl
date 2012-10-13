#version 120
#extension GL_ARB_gpu_shader5 : enable

varying vec2 p;

void main ()
{
    p =  (gl_Vertex).xy;

    gl_Position =   ftransform();
}
