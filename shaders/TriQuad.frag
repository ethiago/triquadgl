#version 120
#extension GL_ARB_gpu_shader_fp64 : enable

varying mat3 Q;
varying vec3 p;

void main ()
{
        gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0) ;

        float f = dot(p,Q*p);

        if(f < 0.0)
            gl_FragColor.r = abs(f);
        else
            gl_FragColor.g = f;

        if(abs(f) < abs(fwidth(f))*0.5)
            gl_FragColor = vec4(1.0) ;
}
