#ifndef CHEBUILDER_H
#define CHEBUILDER_H

#include "compacthalfedge.h"

class CHEBuilder
{

protected:
    CompactHalfEdge m_che;

public:
    CHEBuilder();

    const CompactHalfEdge& che()const;

    virtual void build() = 0;
};

#endif // CHEBUILDER_H
