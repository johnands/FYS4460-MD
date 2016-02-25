#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <string>
#include <vector>
#include "../system.h"

class Potential
{
protected:
    double m_potentialEnergy = 0;
    class System m_system;
public:
    Potential(class System &system);
    virtual ~Potential() {}
    virtual void calculateForces() = 0;
    double potentialEnergy();
    void setPotentialEnergy(double potentialEnergy);
};
#endif
