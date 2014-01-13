#include "fastmarching.h"
#include <QDebug>
#include <qmath.h>

FastMarching::FastMarching(const QImage& src,QObject *parent) :
    QObject(parent), alreadyRun(false)
{
    m_src = src;

    m_values = new float*[m_src.width()];
    for(int i = 0; i < m_src.width(); ++i)
    {
        m_values[i] = new float[m_src.height()];
        for(int j = 0; j < m_src.height(); ++j)
        {
            m_values[i][j] =  -1.0;
        }
    }
}

FastMarching::~FastMarching()
{
    for(int i = 0; i < m_src.width(); ++i)
    {
        delete [] m_values[i];
    }
    delete [] m_values;
}

void FastMarching::allocMap()
{
    m_map    = new bool* [m_src.width()];
    for(int i = 0; i < m_src.width(); ++i)
    {
        m_map   [i] = new bool [m_src.height()];
        for(int j = 0; j < m_src.height(); ++j)
        {
            m_map   [i][j] = false;
        }
    }
}

void FastMarching::freeMap()
{
    for(int i = 0; i < m_src.width(); ++i)
    {
        delete [] m_map[i];
    }
    delete [] m_map;
}

const float ** FastMarching::values()const
{
    return const_cast<const float **>(m_values);
}

float FastMarching::maxValues() const
{
    float max = m_values[0][0];
    for(int i = 0; i < m_src.width(); ++i)
    {
        for(int j = 0; j < m_src.height(); ++j)
        {
            if(m_values[i][j] > max)
                max = m_values[i][j];
        }
    }
    return max;
}

void FastMarching::normalize(float t)
{
    float factor;
    if(t < 0)
        factor = sqrt(2.0)/maxValues();
    else
        factor = 1.0/t;


    for(int i = 0; i < m_src.width(); ++i)
    {
        for(int j = 0; j < m_src.height(); ++j)
        {
            m_values[i][j] *= factor;
        }
    }
}

void FastMarching::run()
{
    if(alreadyRun)
        return;
    alreadyRun = true;

    allocMap();

    for(int j = 0; j < m_src.height(); ++j)
    {
        for(int i = 0; i < m_src.width(); ++i)
        {
            int r = qRed(m_src.pixel(i,j));
            int g = qGreen(m_src.pixel(i,j));
            int b = qBlue(m_src.pixel(i,j));

            if(r < 250 || g < 250 || b < 250)
            {
                fila.push(TexelStruct(i,j,0.0));
                m_map[i][j] = true;
                m_values[i][j] = 0.0;
            }
        }
    }

    flood();

    freeMap();
}

QImage FastMarching::getMapImage()const
{
    QImage result(m_src.size(), QImage::Format_ARGB32);
    float max = maxValues();

    for(int i = 0; i < m_src.width(); ++i)
    {
        for(int j = 0; j < m_src.height(); ++j)
        {
            int v = (int)((m_values[i][j]*255.0)/max + 0.5);
            result.setPixel(i,j,qRgba(v,v,v,255));
        }
    }
    return result;
}

const QImage& FastMarching::getImage()const
{
    return m_src;
}

bool FastMarching::isValid(int i, int j) const
{
    return (i >= 0 && i < m_src.width() && j >=0 && j < m_src.height());
}

void FastMarching::flood()
{
    float d1 = 1.0; float d2 = sqrt(2.0);
    int      i[8] = { -1, -1, -1,  0,  0,  1,  1,  1};
    int      j[8] = { -1,  0,  1, -1,  1, -1,  0,  1};
    float dist[8] = { d2, d1, d2, d1, d1, d2, d1, d2};

    while(!fila.empty())
    {
        TexelStruct curr = fila.top();
        fila.pop();

        int x = curr.i;
        int y = curr.j;

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
                m_values[xi][yj] += dist[k];
                m_values[x][y]    = m_values[xi][yj];
            }
        }
    }
}

float FastMarching::distanceTo(const FastMarching& fm)
{
    float sum = 0.0;
    int cnt = 0;
    for(int i = 0; i < m_src.width(); ++i)
    {
        for(int j = 0; j < m_src.height(); ++j)
        {
            if(m_values[i][j] < 0.1)
            {
                sum += fm.values()[i][j];
                cnt++;
            }
        }
    }
    return sum/cnt;
}
