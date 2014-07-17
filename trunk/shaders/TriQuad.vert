#version 120
#extension GL_ARB_gpu_shader5 : enable

uniform vec2 maior;
uniform vec2 menor;

attribute vec3 abc;
attribute vec3 def;

varying mat3 Q;
varying vec3 P;
varying vec2 texCoord;

varying float k;
varying vec3 corr;

void main ()
{
    Q[0][0]           = abc[0];
    Q[1][0] = Q[0][1] = abc[1];
    Q[2][0] = Q[0][2] = abc[2];

    Q[1][1]           = def[0];
    Q[2][1] = Q[1][2] = def[1];
    Q[2][2]           = def[2];

    P =  (gl_Vertex).xyw;
    k = dot(P,Q*P);
    corr = Q*P;

    gl_Position =   ftransform();
    gl_FrontColor = gl_Color;
    texCoord.x = (P.x - menor.x)/(maior.x - menor.x);
    texCoord.y = (P.y - menor.y)/(maior.y - menor.y);
}
