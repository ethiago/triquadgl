#ifndef FITTING_H
#define FITTING_H

#include <QVector>
#include <QMatrix4x4>
#include <QVector2D>
#include "triquadmesh.h"

QVector<Quadric> fittingGSL(const QMatrix4x4& inv, const QVector<QVector2D>& points);

#endif // FITTING_H
