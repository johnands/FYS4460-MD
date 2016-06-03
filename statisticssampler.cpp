#include "system.h"
#include "statisticssampler.h"
#include "potentials/potential.h"
#include "unitconverter.h"
#include "math/random.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace std;

StatisticsSampler::StatisticsSampler(System *system)
{
    m_system = system;
    m_writeSampleToFile = m_system->getWriteSample();
    m_systemSize = m_system->systemSize().x();
    m_numberOfAtoms = m_system->atoms().size();
}

void StatisticsSampler::saveToFile(int timeStep)
{

    if (m_initialSample) {
        m_outFile.open("argonGasNc10T420Nt501.txt", fstream::out | fstream::trunc);
        m_outFile << setw(12) << "Time step"   << " ";
        m_outFile << setw(12) << "Kinetic"     << " ";
        m_outFile << setw(12) << "Potential"   << " ";
        m_outFile << setw(12) << "Total"       << " ";
        m_outFile << setw(12) << "Temperature" << " ";
        m_outFile << setw(12) << "Pressure"    << endl;

        m_initialSample = false;
        m_outFile.close();
    }

    m_outFile.open("argonGasNc10T420Nt501.txt", fstream::out | fstream::app);

    m_outFile << setw(12) << setprecision(7) << timeStep  << " ";
    m_outFile << setw(12) << setprecision(7) << m_kineticEnergy   << " ";
    m_outFile << setw(12) << setprecision(7) << m_potentialEnergy << " ";
    m_outFile << setw(12) << setprecision(7) << totalEnergy()     << " ";
    m_outFile << setw(12) << setprecision(7) << m_temperature     << " ";
    m_outFile << setw(12) << setprecision(7) << m_pressure        << endl;
    //m_outFile << setw(12) << setprecision(7) << m_meanSquareDisplacement << endl;
    m_outFile.close();
}

void StatisticsSampler::sample(int timeStep)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy();
    samplePotentialEnergy();
    sampleTemperature();
    sampleDensity();
    samplePressure();
    sampleMeanSquareDisplacement();
    if (m_writeSampleToFile) { saveToFile(timeStep); }
}

void StatisticsSampler::sampleKineticEnergy()
{
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
    for(Atom *atom : m_system->atoms()) {
        m_kineticEnergy += 0.5 * atom->mass() * atom->velocity.lengthSquared();
    }
    m_kineticEnergy /= m_numberOfAtoms;
}

void StatisticsSampler::samplePotentialEnergy()
{
    m_potentialEnergy = m_system->potential()->potentialEnergy() / m_numberOfAtoms;
}

void StatisticsSampler::sampleTemperature()
{
    m_temperature = (2*m_kineticEnergy) / 3.0;
}

void StatisticsSampler::samplePressure()
{
    m_pressure = m_system->potential()->pressure() + m_density*m_temperature;
}

void StatisticsSampler::sampleDensity()
{
    m_density = m_numberOfAtoms / (m_systemSize*m_systemSize*m_systemSize);
}

void StatisticsSampler::sampleMeanSquareDisplacement()
{
    m_meanSquareDisplacement = 0;
    for(Atom *atom : m_system->atoms()) {
        m_meanSquareDisplacement += (atom->position - atom->initialPosition()).lengthSquared();
    }
    m_meanSquareDisplacement /= m_numberOfAtoms;
}

void StatisticsSampler::sampleRadialDistribution(int numberOfBins) {
    if (m_firstRadial) {

        m_radialDistFile.open("radialDistribution.txt", fstream::out | fstream::trunc);
        m_radialDistFile << numberOfBins << endl;
        m_radialDistFile << setw(10) << "Bin" << "  ";
        m_radialDistFile << setw(10) << "No. of atoms" << endl;
        m_radialDistFile.close();

        // choose random atom
        m_chosenIndex = floor(Random::nextDouble() * m_numberOfAtoms);
        m_chosenAtom = m_system->atoms()[m_chosenIndex];

        m_radialDistribution.resize(numberOfBins);
        m_bins.resize(numberOfBins+1);

        // make bins
        m_binSize = m_system->systemSizeHalf().x() / numberOfBins;
        for (int i=0; i < numberOfBins; i++) {
            m_bins[i] = i*m_binSize;
            m_radialDistribution[i] = 0;
        }
        m_bins[numberOfBins] = numberOfBins*m_binSize;

        m_firstRadial = false;
    }

    // loop over all pairs of atoms, calculate distance and put in bins
    for (int i=0; i < m_numberOfAtoms; i++) {
        if (i != m_chosenIndex) {
            Atom *atom2 = m_system->atoms()[i];

            // calculate distance
            vec3 dr = m_chosenAtom->position - atom2->position;

            // periodic boundary conditions
            for (int dim=0; dim < 3; dim++) {
                if (dr[dim] > m_system->systemSizeHalf()[dim]) { dr[dim] -= m_system->systemSize()[dim]; }
                else if (dr[dim] < -m_system->systemSizeHalf()[dim]) { dr[dim] += m_system->systemSize()[dim]; }
            }
            double distance = dr.length();

            // fill bins
            for (int bin=0; bin < numberOfBins; bin++) {
                if (distance < m_bins[bin]) {
                    m_radialDistribution[bin] += 1;
                    break;
                }
            }
        }
    }

    // normalize and write to file
    m_radialDistFile.open("radialDistribution.txt", fstream::out | fstream::app);
    for (int bin=0; bin < numberOfBins; bin++) {
        double shellVolume = m_shellVolumeFactor *
        (m_bins[bin+1]*m_bins[bin+1]*m_bins[bin+1] - m_bins[bin]*m_bins[bin]*m_bins[bin]);
        /*if (m_radialDistribution[bin] < shellVolume && m_radialDistribution[bin] > 0) {
            m_radialDistribution[bin] = 1;
        }*/
        m_radialDistribution[bin] /= shellVolume;

        m_radialDistFile << setw(10) << setprecision(4) << m_bins[bin+1] << "  ";
        m_radialDistFile << setw(10) << m_radialDistribution[bin] << endl;
    }

    m_radialDistFile.close();

    // clear radial distribution vector before next computation
    for (int i=0; i < numberOfBins; i++) {
        m_radialDistribution[i] = 0;
    }
}
