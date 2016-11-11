#include "lennardjonescelllist.h"
#include "celllist.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

LennardJonesCellList::LennardJonesCellList(System *system, double sigma, double epsilon,
                                           double rCut, double neighbourCut) :
    Potential (system),
    m_sigma(sigma),
    m_epsilon(epsilon),
    m_sigma6(pow(sigma, 6)),
    m_rCut(rCut),
    m_rCutSquared(rCut*rCut),
    m_neighbourCut(neighbourCut)
{
    m_cellList = new CellList(system, rCut, neighbourCut);
    m_cellList->setupCells();
    m_cellList->updateCells();
    m_updateLists = 0;
    double r2 = 1.0 / m_rCutSquared;
    m_potentialCut = 4*r2*r2*r2*(r2*r2*r2 - 1);
}

void LennardJonesCellList::calculateForces()
{
    double potentialEnergy = 0;
    double pressure = 0;

    // update lists every 5th time step
    if (m_updateLists >= 20) {
        m_cellList->clearCells();
        m_cellList->updateCells();
        m_updateLists = 0;
    }
    m_updateLists++;

    for (int i=0; i < m_system->atoms().size(); i++) {
        Atom *atom1 = m_system->atoms()[i];

        // loop over all atoms in atom1's neighbour list
        for (int j=0; j < m_cellList->getNeighbours()[i].size(); j++) {

            Atom *atom2 = m_cellList->getNeighbours()[i][j];

            // calculate distance vector
            vec3 dr = atom1->position - atom2->position;

            // make sure we're using shortest distance component-wise (periodic boundary conditions)
            for (int dim=0; dim < 3; dim++) {
                if (dr[dim] > m_system->systemSizeHalf()[dim]) { dr[dim] -= m_system->systemSize()[dim]; }
                else if (dr[dim] < -m_system->systemSizeHalf()[dim]) { dr[dim] += m_system->systemSize()[dim]; }
            }

            double dr2 = dr.lengthSquared();
            double r2 = 1.0 / dr2;
            double r6 = r2*r2*r2;

            // check if distance is shorter than cut-off
            // distance could be in range [2.5, 3.0]
            const int cut = (dr2 <= m_rCutSquared);

            // calculate force and potential energy
            vec3 forceOnAtom = 24*r6*(2*r6 - 1)*r2*dr * cut;

            // calculated "shifted and truncated LJ-potential" Vtrunc, which is
            // LJ - LJ_trunc for r < r_cut and 0 otherwise
            potentialEnergy  += ( r6*(r6 - 1) - m_potentialCut ) * cut;
            //cout << forceOnAtom << endl;

            // add contribution to force on atom i and j
            atom1->force += forceOnAtom;
            atom2->force -= forceOnAtom;   // Newton's third law

            if (m_system->getUseExternalForce()) {
                atom1->force[0] += 1.0;
            }

            // dot product of Fij and dr
            pressure += forceOnAtom.dot(dr);
        }
    }

    m_potentialEnergy = 4*potentialEnergy;
    m_pressure = m_inverseVolume*pressure;
}
