#version 120
#extension GL_ARB_gpu_shader5 : enable

attribute vec3 abc;
attribute vec3 def;

varying mat3 Q;
varying vec3 p;

void main ()
{
    Q[0][0]           = abc[0];
    Q[1][0] = Q[0][1] = abc[1];
    Q[2][0] = Q[0][2] = abc[2];

    Q[1][1]           = def[0];
    Q[2][1] = Q[1][2] = def[1];
    Q[2][2]           = def[2];

    p =  (gl_Vertex).xyw;

    gl_Position =   ftransform();
    gl_FrontColor = gl_Color;
}
