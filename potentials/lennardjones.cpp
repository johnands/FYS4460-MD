#include "lennardjones.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

LennardJones::LennardJones(System &system, double sigma, double epsilon) :
    Potential (system),
    m_sigma(sigma),
    m_epsilon(epsilon),
    m_sigma6(pow(sigma, 6))
{

}

void LennardJones::calculateForces()
{
    double potentialEnergy = 0;
    for (int i=0; i < m_system.atoms().size(); i++) {
        vec3 dr, forceOnAtom;
        double dr6, dr2;

        // loop over all other atoms to calculate distances
        for (int j=i+1; j < m_system.atoms().size(); j++) {

            // calculate distance vector
            dr = m_system.atoms()[i]->position - m_system.atoms()[j]->position;

            // make sure we're using shortest distance component-wise (periodic boundary conditions)
            if (dr.x() >  m_system.systemSize().x() / 2.0) { dr[0] -= m_system.systemSize().x(); }
            if (dr.x() < -m_system.systemSize().x() / 2.0) { dr[0] += m_system.systemSize().x(); }

            if (dr.y() >  m_system.systemSize().y() / 2.0) { dr[1] -= m_system.systemSize().y(); }
            if (dr.y() < -m_system.systemSize().y() / 2.0) { dr[1] += m_system.systemSize().y(); }

            if (dr.z() >  m_system.systemSize().z() / 2.0) { dr[2] -= m_system.systemSize().z(); }
            if (dr.z() < -m_system.systemSize().z() / 2.0) { dr[2] += m_system.systemSize().z(); }

            // calculate force
            dr6 = pow(dr.length(), 6);
            dr2 = dr.lengthSquared();
            forceOnAtom = 24*m_epsilon*m_sigma6*(1.0 / dr6)*((2*m_sigma6) / dr6 - 1) * (dr / dr2);

            // add contribution to force on atom i and j
            m_system.atoms()[i]->force += forceOnAtom;
            m_system.atoms()[j]->force -= forceOnAtom;   // Newton's third law

            // calculate potential energy
            potentialEnergy += (m_sigma6 / dr6 - 1) / dr6;
        }
    }


    m_potentialEnergy = potentialEnergy * 4*m_epsilon*m_sigma6;

}
