#ifndef VELOCITYVERLET_H
#define VELOCITYVERLET_H
#include "integrators/integrator.h"

class VelocityVerlet : public Integrator
{
public:
    VelocityVerlet(System *system);
    ~VelocityVerlet() { }
    virtual void integrate(double dt) override;

private:
    bool m_firstStep = true;
};
#endif
