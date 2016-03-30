#ifndef ANDERSEN_H
#define ANDERSEN_H
#include "thermostat.h"
#include "../unitconverter.h"


class Andersen : public Thermostat {
public:
    Andersen(class System *System, double temperature, double dtTau);
    void applyThermostat(double temperature);

private:
    const double m_radiusOfAtom = 2 * UnitConverter::lengthFromSI(71e-12);
};

#endif // ANDERSEN_H
