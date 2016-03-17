#include "thermostat.h"
#include "../system.h"

Thermostat::Thermostat(System &system) {
    m_system = system;
}

void Thermostat::setTemperatureHeatBath(double temperatureHeatBath) {
    m_temperatureHeatBath = temperatureHeatBath;
}

void Thermostat::setDtTau(double dtTau) {
    m_dtTau = dtTau;
}
