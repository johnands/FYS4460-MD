#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H

class System;
class StatisticsSampler
{
private:
    System m_system;
    double m_systemSize;
    double m_numberOfAtoms;

    double m_kineticEnergy = 0;
    double m_potentialEnergy = 0;
    double m_temperature = 0;
    double m_density = 0;
    double m_pressure = 0;
public:
    StatisticsSampler(System &system);
    void saveToFile();
    void sample();
    void sampleKineticEnergy();
    void samplePotentialEnergy();
    void sampleTemperature();
    void samplePressure();
    void sampleDensity();
    double kineticEnergy() { return m_kineticEnergy; }
    double potentialEnergy() { return m_potentialEnergy; }
    double totalEnergy() { return m_kineticEnergy + m_potentialEnergy; }
    double temperature() { return m_temperature; }
    double pressure() { return m_pressure; }
    double density() { return m_density; }
};
#endif
