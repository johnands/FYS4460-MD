#include "lennardjonescelllist.h"
#include "celllist.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

LennardJonesCellList::LennardJonesCellList(System &system, double sigma, double epsilon, double cutOffDistance) :
    Potential (system),
    m_sigma(sigma),
    m_epsilon(epsilon),
    m_sigma6(pow(sigma, 6)),
    m_cutOffDistance(cutOffDistance)
{
    m_cellList = new CellList(system, cutOffDistance);
    m_cellList->setupCells();
}

void LennardJonesCellList::calculateForces()
{
    double potentialEnergy = 0;

    // clear cells before making new ones
    m_cellList->clearCells();
    m_cellList->updateCells();

    for (int i=0; i < m_system.atoms().size(); i++) {
        Atom *atom1 = m_system.atoms()[i];
        vec3 dr, forceOnAtom, cell;
        double dr6, dr2;

        // find out which cell atom i is in
        cell = atom1->cellIndicies();
        cout << cell << endl;

        cout << "neighbour size: " << atom1->neighbourList().size() << endl;

        // loop over all atoms in atom1's neighbour list
        for (int j=0; j < atom1->neighbourList().size(); j++) {
            Atom *atom2 = atom1->neighbourList()[j];

            // calculate distance vector
            dr = atom1->position - atom2->position;

            // make sure we're using shortest distance component-wise (periodic boundary conditions)
            for (int dim=0; dim < 3; dim++) {
                if (dr[dim] > m_system.systemSizeHalf()[dim]) { dr[dim] -= m_system.systemSize()[dim]; }
                else if (dr[dim] < -m_system.systemSizeHalf()[dim]) { dr[dim] += m_system.systemSize()[dim]; }
            }

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
