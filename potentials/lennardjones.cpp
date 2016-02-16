#include "lennardjones.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon),
    m_sigma6(pow(sigma, 6))
{

}

void LennardJones::calculateForces(System *system)
{
    double potentialEnergy = 0;
    for (int i=0; i < system->atoms().size(); i++) {
        vec3 dr, forceOnAtom;
        double dr6, dr2;

        // loop over all other atoms to calculate distances
        for (int j=i+1; j < system->atoms().size(); j++) {

            // calculate distance vector
            dr = system->atoms()[i]->position - system->atoms()[j]->position;

            // make sure we're using shortest distance component-wise (periodic boundary conditions)
            if (dr.x() >  system->systemSize().x() / 2) { dr[0] -= system->systemSize().x(); }
            if (dr.x() < -system->systemSize().x() / 2) { dr[0] += system->systemSize().x(); }

            if (dr.y() >  system->systemSize().y() / 2) { dr[1] -= system->systemSize().y(); }
            if (dr.y() < -system->systemSize().y() / 2) { dr[1] += system->systemSize().y(); }

            if (dr.z() >  system->systemSize().z() / 2) { dr[2] -= system->systemSize().z(); }
            if (dr.z() < -system->systemSize().z() / 2) { dr[2] += system->systemSize().z(); }

            // calculate force
            dr6 = pow(dr.length(), 6);
            dr2 = dr.lengthSquared();
            forceOnAtom = -24*m_epsilon*m_sigma6*(1 / dr6)*((2*m_sigma6) / dr6 - 1) * (dr / dr2);

            // add contribution to force on atom i and j
            system->atoms()[i]->force += forceOnAtom;
            system->atoms()[j]->force -= forceOnAtom;   // Newton's third law

            // calculate potential energy
            potentialEnergy += (m_sigma6 / dr6 - 1) / dr6;
        }
    }


    m_potentialEnergy = potentialEnergy * 4*m_epsilon*m_sigma6;

}
