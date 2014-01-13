#ifndef ERRORCOMPUTING_H
#define ERRORCOMPUTING_H

#include <QObject>
#include <QImage>
#include "fastmarching.h"

class ErrorComputing : public QObject
{
    Q_OBJECT

    FastMarching *m_input;
    FastMarching *m_triquad;

public:
    explicit ErrorComputing(const QImage& input, const QImage& triquad);
    ~ErrorComputing();

    float getNormalizedFittingError();
    float getNormalizedExtraConponentsError();
    
signals:
    
public slots:
    
};

#endif // ERRORCOMPUTING_H
