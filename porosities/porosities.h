#ifndef POROSITIES_H
#define POROSITIES_H
#include "../system.h"

class Porosities
{
public:
    Porosities(class System &system);
    virtual void makePores() = 0;

protected:
    class System m_system;
};

#endif // POROSITIES_H
