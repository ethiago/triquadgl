#extension GL_ARB_gpu_shader5:enable

varying vec3 p;

vec3 v1 = vec3(-1.0, -0.5, 1.0);
vec3 v2 = vec3( 1.0, -0.5, 1.0);
vec3 v3 = vec3( 0.0,  0.5, 1.0);

mat3 q1 = mat3( 4.0, 0.0, 0.0,   0.0, 0.0, 0.5,   0.0, 0.5, 0.0);
mat3 q2 = mat3(-4.0, 0.0, 0.0,   0.0, 0.0, 0.5,   0.0, 0.5, 0.0);
mat3 q3 = mat3( 1.0, 0.0, 0.0,   0.0, 1.0, 0.0,   0.0, 0.0, 1.0);

float minimum_distance(vec2 v, vec2 w, vec2 p) 
{
  float L2 = dot(v-w,v-w);

  float t = dot(p - v, w - v) / L2;
  if (t < 0.0) return distance(p, v);       // Beyond the 'v' end of the segment
  else if (t > 1.0) return distance(p, w);  // Beyond the 'w' end of the segment
  vec2 projection = v + t * (w - v);  // Projection falls on the segment
  return distance(p, projection);
}
void main()
{
	mat3 W = mat3(v1,v2,v3);
	W = inverse(W);
	vec3 b = W*p;
	mat3 q = b[0]*q1 + b[1]*q2 + b[2]*q3;

	float f =  dot(p, q*p);
	float d1 = minimum_distance(v1.xy, v2.xy, p.xy);
	float d2 = minimum_distance(v2.xy, v3.xy, p.xy);
	float d3 = minimum_distance(v3.xy, v1.xy, p.xy);

	gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);

	if(abs(f) < fwidth(f)*1.0)
		gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);
	else if ( abs(d1) < fwidth(d1) ||  abs(d2) < fwidth(d2) || abs(d3) < fwidth(d3) )
		gl_FragColor = vec4(0.5, 0.5, 0.5, 1.0);
	else
		discard;
}