uniform sampler2D sampler2d0;
uniform vec2 maior;
uniform vec2 menor;

uniform int showScalar;

varying vec2 texCoord;

varying mat3 Qv;
varying vec3 Po;

varying mat3 Qx;
varying mat3 Qy;



//Dados de Entrada--------------------------------------------------------------------------------
const ivec2 dimensaoRuido = ivec2(1200,1200);
const int L = int(dimensaoRuido.x/20);
float dx = 1.0/float(dimensaoRuido.x);//tamanho do pixel em x
float dy = 1.0/float(dimensaoRuido.y);//tamanho do pixel em y	


vec2 mapToDomain(vec2 texcoord)
{
	return vec2(texcoord.x*(maior.x-menor.x)+menor.x, texcoord.y*(maior.y-menor.y)+menor.y);
}

//Funcoes-----------------------------------------------------------------------------------------
//Defenido o Campo Vetorial
vec2 campoVetorial(vec2 ponto)
{
       vec3 T =2.0*(Qv*vec3(ponto,1.0));
        return vec2( T.x + dot(Po,Qx*Po), T.y + dot(Po,Qy*Po));
}

//Computando a Covolucao ao longo da Curva Integral reverente ao pixel ponto=(x,y)=texCoord
vec3 convolucaoCurvaIntegral (vec2 ponto)
{
	vec2 p=ponto;
	vec2 V=normalize(campoVetorial(mapToDomain(p)));
	vec3 somaCorCurvaIntegral=vec3(0.0);
	float ds=min(dx,dy);
	somaCorCurvaIntegral+=texture2D(sampler2d0,p).rgb;
	for(int s=0;s<L;s=s+1)//calculando os pontos da linha de fluxo no sentido positivo
	{
		p+=ds*V;
		somaCorCurvaIntegral+=texture2D(sampler2d0,p).rgb;
		V=normalize(campoVetorial(mapToDomain(p)));
	}
	p=ponto;
	V=normalize(campoVetorial(mapToDomain(p)));
	for(int s=0;s>-L;s=s-1)//calculando os pontos da linha de fluxo no sentido negativo
	{
		p-=ds*V;
		somaCorCurvaIntegral+=texture2D(sampler2d0,p).rgb*float(-s);
		V=normalize(campoVetorial(mapToDomain(p)));
	}
	somaCorCurvaIntegral/=float(0.5*(1.0+float(L))*float(L) + float(L) +1.0);

	return somaCorCurvaIntegral;
}



//Main--------------------------------------------------------------------------------------------
void main()
{
    float prop = 0.7;
        float f = dot(Po,Qv*Po);
	if(abs(f) < abs(fwidth(f))*0.5)
	{
                gl_FragColor = vec4(1.0) ;
	}
	else
	{
            vec3 cor = convolucaoCurvaIntegral(texCoord);
            if(showScalar == 0)
            {
                gl_FragColor.rgb = cor;
            }else
            {
                if(f < 0.0)
                {
                    f = min(-f, 1.0);
                    gl_FragColor.rgb = cor*prop + vec3( f,0.0,0.0)*(1.0-prop);
                }else
                {
                    f = min(f, 1.0);
                    gl_FragColor.rgb = cor*prop + vec3(0.0,f,0.0)*(1.0-prop);
                }
            }
            gl_FragColor.a = 1.0;
	}
}
