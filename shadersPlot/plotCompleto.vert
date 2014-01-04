
varying vec3 p;

void main()
{
	vec4 v = gl_Vertex;
	v.xyz *= 2.0;

	p = v.xyw;
	gl_Position = gl_ModelViewProjectionMatrix * v;
}