#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H
#include <fstream>
#include <vector>

class System;
class StatisticsSampler
{
private:
    System *m_system = nullptr;
    double m_systemSize;
    double m_numberOfAtoms;
    std::fstream m_outFile;
    std::fstream m_radialDistFile;
    bool m_initialSample = true;

    std::vector<double> m_radialDistribution;
    std::vector<double> m_bins;
    double m_binSize;
    bool m_firstRadial = true;
    int m_chosenIndex;
    class Atom *m_chosenAtom = nullptr;
    const double m_shellVolumeFactor = 4.0*3.14159265359/3.0;
    bool m_writeSampleToFile;

    double m_kineticEnergy = 0;
    double m_kineticEnergyPerAtom = 0;
    double m_potentialEnergy = 0;
    double m_potentialEnergyPerAtom = 0;
    double m_temperature = 0;
    double m_density = 0;
    double m_pressure = 0;
    double m_meanSquareDisplacement = 0;
    double m_flowVelocity = 0;

public:
    StatisticsSampler(System *system);
    void saveToFile(int timeStep);
    void sample(int timeStep);
    void sampleKineticEnergy();
    void samplePotentialEnergy();
    void sampleTemperature();
    void samplePressure();
    void sampleDensity();
    void sampleMeanSquareDisplacement();
    void sampleRadialDistribution(int numberOfBins);
    void sampleFlowVelocity();

    double kineticEnergy() { return m_kineticEnergy; }
    double kineticEnergyPerAtom() { return m_kineticEnergyPerAtom; }
    double potentialEnergy() { return m_potentialEnergy; }
    double potentialEnergyPerAtom() { return m_potentialEnergyPerAtom; }
    double totalEnergy() { return m_kineticEnergy + m_potentialEnergy; }
    double totalEnergyPerAtom() { return m_kineticEnergyPerAtom + m_potentialEnergyPerAtom; }
    double temperature() { return m_temperature; }
    double pressure() { return m_pressure; }
    double density() { return m_density; }
    double MeanSquareDisplacement() { return m_meanSquareDisplacement; }
    double flowVelocity() { return m_flowVelocity; }
};
#endif
