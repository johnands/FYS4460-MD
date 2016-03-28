#include "thermostat.h"
#include "../system.h"
#include <iostream>

Thermostat::Thermostat(System &system, double temperatureHeatBath, double tau) {
    m_system = system;
    m_temperatureHeatBath = temperatureHeatBath;
    m_dtTau = m_system.getTimeStep()/tau;
}

void Thermostat::setTemperatureHeatBath(double temperatureHeatBath) {
    m_temperatureHeatBath = temperatureHeatBath;
}

void Thermostat::setDtTau(double dtTau) {
    m_dtTau = dtTau;
}
