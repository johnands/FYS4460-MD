#ifndef BERENDSEN_H
#define BERENDSEN_H
#include "thermostat.h"


class Berendsen : public Thermostat {
public:
    Berendsen(class System *System, double temperature, double dtTau);
    void applyThermostat(double temperature);
};

#endif // BERENDSEN_H
