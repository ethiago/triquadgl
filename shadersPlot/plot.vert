varying float t;
varying vec3 p;

void main()
{
	t = (gl_Vertex.x + 1.0)/2.0;
	p = gl_Vertex.xyw;
	gl_Position = ftransform();
}
