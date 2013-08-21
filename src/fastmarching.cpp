#include "fastmarching.h"
#include <QDebug>
#include <qmath.h>

FastMarching::FastMarching(const QImage& src,QObject *parent) :
    QObject(parent)
{
    m_src = src;
}

void FastMarching::setSource(const QImage& src)
{
    m_src = src;
}

void FastMarching::allocMats()
{
    map    = new bool* [m_src.width()];
    values = new float*[m_src.width()];
    for(int i = 0; i < m_src.width(); ++i)
    {
        map   [i] = new bool [m_src.height()];
        values[i] = new float[m_src.height()];
        for(int j = 0; j < m_src.height(); ++j)
        {
            map   [i][j] = false;
            values[i][j] =  -1.0;
        }
    }
}

void FastMarching::freeMats()
{
    for(int i = 0; i < m_src.width(); ++i)
    {
        delete [] map[i];
        delete [] values[i];
    }
    delete [] map;
    delete [] values;
}

float FastMarching::maxValues()
{
    float max = values[0][0];
    for(int i = 0; i < m_src.width(); ++i)
    {
        for(int j = 0; j < m_src.height(); ++j)
        {
            if(values[i][j] > max)
                max = values[i][j];
        }
    }
    return max;
}

QImage FastMarching::run()
{
    QImage result(m_src.size(), m_src.format());

    allocMats();

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
                map[i][j] = true;
                values[i][j] = 0.0;
            }
        }
    }

    flood();


    float max = maxValues();

    for(int i = 0; i < m_src.width(); ++i)
    {
        for(int j = 0; j < m_src.height(); ++j)
        {
            int v = (int)((values[i][j]*255.0)/max + 0.5);
            result.setPixel(i,j,qRgba(v,v,v,255));
        }
    }


    freeMats();

    return result;
}

bool FastMarching::isValid(int i, int j)
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

        //qDebug() <<QPoint(x,y);
        for(int k = 0; k < 8; ++k)
        {
            int xi = x+i[k];
            int yj = y+j[k];
            //qDebug() << " " <<QPoint(x+i,y+j);
            if(!isValid(xi,yj))
                continue;

            if(!map[xi][yj])
            {
                values[xi][yj] = curr.value + dist[k];
                map[xi][yj] = true;
                fila.push(TexelStruct(xi,yj,values[xi][yj]));
            }else if(curr.value + dist[k] < values[xi][yj])
            {
                values[xi][yj] = curr.value + dist[k];
            }else if(values[xi][yj] + dist[k]  <  curr.value)
            {
                values[x][y]    = values[xi][yj] + dist[k];
                values[xi][yj] += dist[k];
            }
        }
    }
}
