#include "fakefm.h"
#include <QDebug>
#include <qmath.h>
#include <limits>

#define INF 1e+08

FakeFM::FakeFM(QObject *parent) :
    QObject(parent), alreadyRun(false)
{
    w = 8, h = 8;

    m_values = new float*[w];
    for(int i = 0; i < w; ++i)
    {
        m_values[i] = new float[h];
        for(int j = 0; j < h; ++j)
        {
            m_values[i][j] =  INF;
        }
    }
}

FakeFM::~FakeFM()
{
    for(int i = 0; i < w; ++i)
    {
        delete [] m_values[i];
    }
    delete [] m_values;
}

void FakeFM::allocMap()
{
    m_map    = new bool* [w];
    for(int i = 0; i < w; ++i)
    {
        m_map   [i] = new bool [h];
        for(int j = 0; j < h; ++j)
        {
            m_map   [i][j] = false;
        }
    }
}

void FakeFM::freeMap()
{
    for(int i = 0; i < w; ++i)
    {
        delete [] m_map[i];
    }
    delete [] m_map;
}

void FakeFM::run()
{
    if(alreadyRun)
        return;
    alreadyRun = true;

    allocMap();

    int qtd = 14;
    int I[] = {0,0,1,1,2,3,4,4,5,5,6,6,6,6};
    int J[] = {0,1,1,2,2,2,2,3,3,4,4,5,6,7};

    for(int k = 0; k < qtd; ++k)
    {
        fila.push(TexelStruct(I[k],J[k],0.0));
        m_map[I[k]][J[k]] = true;
        m_values[I[k]][J[k]] = 0.0;
    }

    flood();

    freeMap();
}

bool FakeFM::isValid(int i, int j) const
{
    return (i >= 0 && i < w && j >=0 && j < h);
}

void FakeFM::flood()
{
    float d1 = 1.0; float d2 = sqrt(2.0);
    int      i[8] = { -1, -1, -1,  0,  0,  1,  1,  1};
    int      j[8] = { -1,  0,  1, -1,  1, -1,  0,  1};
    float dist[8] = { d2, d1, d2, d1, d1, d2, d1, d2};


    QDebug deb = qDebug();
    for(int j = 0; j < h; ++j)
    {
        for(int i = 0; i < w; ++i)
        {
            deb << m_values[i][j] << "\t";
        }
        deb << "\n";
    }
    deb << "\n";

    while(!fila.empty())
    {
        TexelStruct curr = fila.top();
        fila.pop();

        int x = curr.i;
        int y = curr.j;
        float a = 0.0;

        for(int k = 0; k < 8; ++k)
        {
            int xi = x+i[k];
            int yj = y+j[k];
            if(!isValid(xi,yj))
                continue;

            if(!m_map[xi][yj])
            {
                m_values[xi][yj] = curr.value + dist[k];
                m_map[xi][yj] = true;
                fila.push(TexelStruct(xi,yj,m_values[xi][yj]));
            }else if(curr.value + dist[k] < m_values[xi][yj])
            {
                m_values[xi][yj] = curr.value + dist[k];
            }else if(m_values[xi][yj] + dist[k]  <  curr.value)
            {
                curr.value = m_values[xi][yj] + dist[k];
                m_values[x][y]    = curr.value;
                m_values[xi][yj] += dist[k];
            }
        }

    }

    for(int j = 0; j < h; ++j)
    {
        for(int i = 0; i < w; ++i)
        {
            deb << m_values[i][j] << "\t";
        }
        deb << "\n";
    }
}

