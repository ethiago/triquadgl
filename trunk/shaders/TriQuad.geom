#version 150
#extension GL_EXT_geometry_shader4 : enable

in mat3 Q[3];
in vec3 P[3];
in vec2 texCoordIn[3];

out mat3 Qv;
out vec3 Po;
out vec2 texCoord;

out mat3 Qx;
out mat3 Qy;

void main()
{
    mat3 W;
    for(int i = 0; i < 3; ++i)
    {
        W[i] = P[i];
    }

    W = inverse(W);
    Qx = Q[0]*W[0][0] + Q[1]*W[0][1] + Q[2]*W[0][2];
    Qy = Q[0]*W[1][0] + Q[1]*W[1][1] + Q[2]*W[1][2];

    for(int i = 0; i < 3; ++i)
    {
        gl_FrontColor = gl_FrontColorIn[i];
        gl_Position = gl_PositionIn[i];
        Qv = Q[i];
        Po = P[i];
        texCoord = texCoordIn[i];
        EmitVertex();
    }
    EndPrimitive();
}
