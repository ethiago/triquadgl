varying float t;
varying vec3 p;

mat3 Q1 = mat3(0.0, 0.0, 0.5,  0.0, 2.0, 0.0,  0.5, 0.0, 0.0);
mat3 Q2 = mat3(0.0, 0.0, 0.5,  0.0,-3.0, 0.0,  0.5, 0.0, 0.0);



void main()
{
	float s = t;

	mat3 Q = (1.0 - s)*Q1 + s*Q2;

	float f = dot(p, Q*p);

	if(abs(f) > fwidth(f) )
		discard;

	gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);
}