#ifndef FITTING_H
#define FITTING_H

#include <QVector>
#include <QMatrix4x4>
#include <QVector2D>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "triquadmesh.h"
#include "quadric2d.h"

QVector<Quadric2D> fittingGSL(const QMatrix4x4& inv, const QVector<QVector2D>& points);
QVector<Quadric2D> fittingGLOBAL(gsl_matrix * A, gsl_vector * B, float f);
QVector<Quadric2D> fittingGLOBAL_flivre(gsl_matrix * A, gsl_vector * B);

#endif // FITTING_H
