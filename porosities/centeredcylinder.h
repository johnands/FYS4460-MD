#ifndef CENTEREDCYLINDER_H
#define CENTEREDCYLINDER_H
#include "porosities.h"
#include "../unitconverter.h"

class CenteredCylinder : public Porosities {
public:
    CenteredCylinder(class System *system, bool usePores);
    void makePores();

private:
    const double m_poreRadius = 2e-9/1e-10;

};

#endif // CENTEREDCYLINDER_H
