#include "errorcomputing.h"

ErrorComputing::ErrorComputing(const QImage& input, const QImage& triquad) :
    QObject(NULL)
{
    int minX = qMax(input.width(), triquad.width());
    int minY = qMax(input.height(), triquad.height());

    int maxX = 0, maxY = 0;

    for(int i = 0 ; i < input.width(); ++i)
    {
        for(int j = 0; j < input.height(); ++j)
        {

            int r = qRed(input.pixel(i,j));
            int g = qGreen(input.pixel(i,j));
            int b = qBlue(input.pixel(i,j));

            if(r < 250 && g < 250 && b < 250)
            {
                if(i < minX)
                    minX = i;
                if(j < minY)
                    minY = j;
                if(i > maxX)
                    maxX = i;
                if(j > maxY)
                    maxY = j;
            }

        }
    }

    for(int i = 0 ; i < triquad.width(); ++i)
    {
        for(int j = 0; j < triquad.height(); ++j)
        {

            int r = qRed(triquad.pixel(i,j));
            int g = qGreen(triquad.pixel(i,j));
            int b = qBlue(triquad.pixel(i,j));

            if(r < 250 && g < 250 && b < 250)
            {
                if(i < minX)
                    minX = i;
                if(j < minY)
                    minY = j;
                if(i > maxX)
                    maxX = i;
                if(j > maxY)
                    maxY = j;
            }

        }
    }

    m_input   = new FastMarching(  input.copy(minX, minY, maxX-minX, maxY-minY));
    m_triquad = new FastMarching(triquad.copy(minX, minY, maxX-minX, maxY-minY));

    m_input->getImage().save("input.png");
    m_triquad->getImage().save("triquad.png");

    m_input->run();
    m_triquad->run();

    m_input->getMapImage().save("inputMap.png");
    m_triquad->getMapImage().save("triquadMap.png");
}

ErrorComputing::~ErrorComputing()
{
    delete m_input;
    delete m_triquad;
}

float ErrorComputing::getNormalizedFittingError()
{
    int norm = qMax(m_input->getImage().width(), m_input->getImage().height());
    return m_input->distanceTo(*m_triquad)/norm;
}

float ErrorComputing::getNormalizedExtraConponentsError()
{
    int norm = qMax(m_input->getImage().width(), m_input->getImage().height());
    return m_triquad->distanceTo(*m_input)/norm;
}
