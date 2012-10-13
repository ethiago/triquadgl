#version 120
#extension GL_ARB_gpu_shader_fp64 : enable

uniform float consts[10];

varying vec2 p;

float f(vec2 v)
{
    return  consts[0]*v.x*v.x*v.x + consts[1]*v.x*v.x*v.y + consts[2]*v.x*v.y*v.y + consts[3]*v.y*v.y*v.y +
            consts[4]*v.x*v.x     + consts[5]*v.x*v.y     + consts[6]*v.y*v.y     +
            consts[7]*v.x         + consts[8]*v.y         +
            consts[9];
}

void main ()
{
        float val = f(p);

        gl_FragColor = vec4(1.0) ;

        if(abs(val) > abs(fwidth(val))/2.0)
            discard;
}
