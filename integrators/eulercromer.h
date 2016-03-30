#ifndef EULERCROMER_H
#define EULERCROMER_H
#include "integrators/integrator.h"

class EulerCromer : public Integrator
{
public:
    EulerCromer(System *system);
    ~EulerCromer() {}
    virtual void integrate(double dt) override;
};

#endif
