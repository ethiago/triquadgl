#version 120
#extension GL_ARB_gpu_shader_fp64 : enable

varying mat3 Q;
varying vec3 p;

void main ()
{
        gl_FragColor = vec4(1.0) ;

        if(abs(dot(p,Q*p)) > abs(fwidth(dot(p,Q*p)))/2.0)
            discard;
}
