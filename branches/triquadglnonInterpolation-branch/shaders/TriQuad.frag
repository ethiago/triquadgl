#version 120
#extension GL_ARB_gpu_shader_fp64 : enable

uniform mat3 Q0;
uniform mat3 Q1;
uniform mat3 Q2;

varying vec3 p;

void main ()
{
        gl_FragColor = vec4(1.0) ;

        mat3 Q = Q0*p.x + Q1*p.y + Q2*p.z;

        if(abs(dot(p,Q*p)) > abs(fwidth(dot(p,Q*p)))*0.5)
            discard;
}
