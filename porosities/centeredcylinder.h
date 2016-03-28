#ifndef CENTEREDCYLINDER_H
#define CENTEREDCYLINDER_H
#include "porosities.h"
#include "../unitconverter.h"

class CenteredCylinder : public Porosities {
public:
    CenteredCylinder(class System &system);
    void makePores();

private:
    const double m_poreRadius = UnitConverter::lengthFromSI(2e-9);
};

#endif // CENTEREDCYLINDER_H
