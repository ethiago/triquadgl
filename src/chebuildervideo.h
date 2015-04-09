#ifndef CHEBUILDERVIDEO_H
#define CHEBUILDERVIDEO_H

#include "chebuilder.h"

class CHEBuilderVideo : public CHEBuilder
{
    int m_type;
public:
    CHEBuilderVideo(int);

    virtual void build();
};

#endif // CHEBUILDERVIDEO_H
