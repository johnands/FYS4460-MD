#include "lennardjones.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

LennardJones::LennardJones(System *system, double sigma, double epsilon) :
    Potential (system),
    m_sigma(sigma),
    m_epsilon(epsilon),
    m_sigma6(pow(sigma, 6))
{

}

void LennardJones::calculateForces()
{
    double potentialEnergy = 0;
    double pressure = 0;
    for (int i=0; i < m_system->atoms().size(); i++) {
        vec3 dr, forceOnAtom;
        double dr6, dr2;

        Atom *atom1 = m_system->atoms()[i];

        // loop over all other atoms to calculate distances
        for (int j=i+1; j < m_system->atoms().size(); j++) {

            Atom *atom2 = m_system->atoms()[j];

            // calculate distance vector
            dr = atom1->position - atom2->position;

            // make sure we're using shortest distance component-wise (periodic boundary conditions)
            if (dr.x() >  m_system->systemSizeHalf().x()) { dr[0] -= m_system->systemSize().x(); }
            if (dr.x() < -m_system->systemSizeHalf().x()) { dr[0] += m_system->systemSize().x(); }

            if (dr.y() >  m_system->systemSizeHalf().y()) { dr[1] -= m_system->systemSize().y(); }
            if (dr.y() < -m_system->systemSizeHalf().y()) { dr[1] += m_system->systemSize().y(); }

            if (dr.z() >  m_system->systemSizeHalf().z()) { dr[2] -= m_system->systemSize().z(); }
            if (dr.z() < -m_system->systemSizeHalf().z()) { dr[2] += m_system->systemSize().z(); }

            // calculate force
            dr6 = 1.0 / pow(dr.length(), 6);
            dr2 = 1.0 / dr.lengthSquared();
            forceOnAtom = 24*dr6*(2*dr6 - 1) * dr2*dr;

            // add contribution to force on atom i and j
            atom1->force += forceOnAtom;
            atom2->force -= forceOnAtom;   // Newton's third law

            // dot product of Fij and dr
            pressure += forceOnAtom.dot(dr);

            // calculate potential energy
            potentialEnergy += (dr6 - 1)*dr6;
        }
    }


    m_potentialEnergy = potentialEnergy * 4;
    m_pressure = m_inverseVolume*pressure;
}
