#ifndef POROSITIES_H
#define POROSITIES_H
#include "../system.h"

class Porosities
{
public:
    Porosities(class System &system);
    void computePorosity();
    void halfDensity();
    virtual void makePores() = 0;

protected:
    class System m_system;
    double m_porosity = 0;
    int m_numberOfMovingAtoms = 0;
};

#endif // POROSITIES_H
