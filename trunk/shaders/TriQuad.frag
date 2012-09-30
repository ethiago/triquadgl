#version 120
#extension GL_ARB_gpu_shader_fp64 : enable

varying mat3 Q;
varying vec3 p;

const float err = 0.01;

void main ()
{
        gl_FragColor = vec4(1.0);

        if(abs(dot(p,Q*p)) > err)
            discard;
}
