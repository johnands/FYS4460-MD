#include "potential.h"

Potential::Potential(System &system) {
    m_system = system;
    m_inverseVolume = 1.0 / ( 3 * m_system.systemSize().x()*m_system.systemSize().y()*m_system.systemSize().z() );
}

double Potential::potentialEnergy() {
    return m_potentialEnergy;
}

double Potential::pressure() {
    return m_pressure;
}

void Potential::setPotentialEnergy(double potentialEnergy) {
    m_potentialEnergy = potentialEnergy;
}

void Potential::setPressure(double pressure) {
    m_pressure = pressure;
}
