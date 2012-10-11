#version 120
#extension GL_ARB_gpu_shader5 : enable

varying vec3 p;

void main ()
{
    p =  (gl_Vertex).xyw;

    gl_Position =   ftransform();
    gl_FrontColor = gl_Color;
}
