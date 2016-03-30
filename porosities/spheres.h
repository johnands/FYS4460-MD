#ifndef SPHERES_H
#define SPHERES_H
#include "porosities.h"
#include "../unitconverter.h"

class Spheres : public Porosities {
public:
    Spheres(System *system);
    void makePores();

private:
    int m_numberOfPores = 20;
    const double m_radius0 = UnitConverter::lengthFromSI(2e-9);
    const double m_radius1 = UnitConverter::lengthFromSI(3e-9);
};

#endif // SPHERES_H
