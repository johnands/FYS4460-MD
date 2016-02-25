#include "potential.h"

Potential::Potential(System &system)
{
    m_system = system;
}

double Potential::potentialEnergy()
{
    return m_potentialEnergy;
}

void Potential::setPotentialEnergy(double potentialEnergy)
{
    m_potentialEnergy = potentialEnergy;
}
