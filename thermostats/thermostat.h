#ifndef THERMOSTAT_H
#define THERMOSTAT_H
#include "../system.h"

class Thermostat
{
public:
    Thermostat(class System &system);
    virtual void applyThermostat(double temperature) = 0;
    void setTemperatureHeatBath(double temperatureHeatBath);
    void setDtTau(double dtTau);

protected:
    class System m_system;
    double m_temperatureHeatBath = 0;
    double m_dtTau = 0;
};

#endif // THERMOSTAT_H
