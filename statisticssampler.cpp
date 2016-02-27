#include "system.h"
#include "statisticssampler.h"
#include "potentials/potential.h"
#include "unitconverter.h"

StatisticsSampler::StatisticsSampler(System &system)
{
    m_system = system;
    m_systemSize = m_system.systemSize().x();
    m_numberOfAtoms = m_system.atoms().size();
}

void StatisticsSampler::saveToFile()
{
    // Save the statistical properties for each timestep for plotting etc.
}

void StatisticsSampler::sample()
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy();
    samplePotentialEnergy();
    sampleTemperature();
    sampleDensity();
    samplePressure();
    saveToFile();
}

void StatisticsSampler::sampleKineticEnergy()
{
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
    for(Atom *atom : m_system.atoms()) {
        m_kineticEnergy += 0.5 * atom->mass() * atom->velocity.lengthSquared();
    }
    m_kineticEnergy /= m_numberOfAtoms;
}

void StatisticsSampler::samplePotentialEnergy()
{
    m_potentialEnergy = m_system.potential()->potentialEnergy() / m_numberOfAtoms;
}

void StatisticsSampler::sampleTemperature()
{
    m_temperature = UnitConverter::temperatureToSI( (2*m_kineticEnergy) / 3.0 );
}

void StatisticsSampler::samplePressure()
{
    m_pressure = m_system.potential()->pressure() + m_density*m_temperature;
}

void StatisticsSampler::sampleDensity()
{
    m_density = m_numberOfAtoms / (m_systemSize*m_systemSize*m_systemSize);
}
