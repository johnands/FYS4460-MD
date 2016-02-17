#include "lennardjonescelllist.h"
#include "celllist.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

LennardJonesCellList::LennardJonesCellList(double sigma, double epsilon, double cutOffDistance) :
    m_sigma(sigma),
    m_epsilon(epsilon),
    m_sigma6(pow(sigma, 6)),
    m_cutOffDistance(cutOffDistance)
{
    m_cellList = new CellList(cutOffDistance);
}

void LennardJonesCellList::calculateForces(System *system)
{
    double potentialEnergy = 0;

    // clear cells before making new ones
    m_cellList->clearCells();
    m_cellList->setup(system);

    for (int i=0; i < system->atoms().size(); i++) {
        vec3 dr, forceOnAtom, cell;
        double dr6, dr2;

        // find out which cell atom i is in
        cell = system->atoms()[i]->cellIndicies;
        cout << cell;

        // loop over all other atoms to calculate distances
        for (int j=i+1; j < system->atoms().size(); j++) {

            // calculate distance vector
            dr = system->atoms()[i]->position - system->atoms()[j]->position;

            // make sure we're using shortest distance component-wise (periodic boundary conditions)
            if (dr.x() >  system->systemSize().x() / 2.0) { dr[0] -= system->systemSize().x(); }
            if (dr.x() < -system->systemSize().x() / 2.0) { dr[0] += system->systemSize().x(); }

            if (dr.y() >  system->systemSize().y() / 2.0) { dr[1] -= system->systemSize().y(); }
            if (dr.y() < -system->systemSize().y() / 2.0) { dr[1] += system->systemSize().y(); }

            if (dr.z() >  system->systemSize().z() / 2.0) { dr[2] -= system->systemSize().z(); }
            if (dr.z() < -system->systemSize().z() / 2.0) { dr[2] += system->systemSize().z(); }

            // calculate force
            dr6 = pow(dr.length(), 6);
            dr2 = dr.lengthSquared();
            forceOnAtom = 24*m_epsilon*m_sigma6*(1.0 / dr6)*((2*m_sigma6) / dr6 - 1) * (dr / dr2);

            // add contribution to force on atom i and j
            system->atoms()[i]->force += forceOnAtom;
            system->atoms()[j]->force -= forceOnAtom;   // Newton's third law

            // calculate potential energy
            potentialEnergy += (m_sigma6 / dr6 - 1) / dr6;
        }
    }


    m_potentialEnergy = potentialEnergy * 4*m_epsilon*m_sigma6;

}
