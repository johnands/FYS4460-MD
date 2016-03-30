#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <string>
#include <vector>
#include "../system.h"

class Potential
{
protected:
    double m_potentialEnergy = 0;
    double m_pressure = 0;
    double m_inverseVolume;
    class System *m_system = nullptr;
public:
    Potential(class System *system);
    virtual ~Potential() {}
    virtual void calculateForces() = 0;
    double potentialEnergy();
    double pressure();
    void setPotentialEnergy(double potentialEnergy);
    void setPressure(double pressure);
};
#endif
