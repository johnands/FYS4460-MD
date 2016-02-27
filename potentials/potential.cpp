#include "potential.h"

Potential::Potential(System &system) {
    m_system = system;
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
