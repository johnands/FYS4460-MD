#ifndef INTEGRATOR_H
#define INTEGRATOR_H

class System;
class Integrator
{
public:
    Integrator(System *system);
    virtual ~Integrator() { }
    virtual void integrate(double dt) = 0;

protected:
    System *m_system = nullptr;
};

#endif
