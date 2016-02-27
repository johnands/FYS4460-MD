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
    m_cellList->updateCells();
    m_updateLists = 0;
    m_inverseVolume = 1.0 / ( 3 * m_system.systemSize().x()*m_system.systemSize().y()*m_system.systemSize().z() );
}

void LennardJonesCellList::calculateForces()
{
    double potentialEnergy = 0;
    double pressure = 0;

    m_cellList->clearCells();
    m_cellList->updateCells();

    // update lists every 10th time step


    for (int i=0; i < m_system.atoms().size(); i++) {
        Atom *atom1 = m_system.atoms()[i];
        vec3 dr, forceOnAtom, cell;
        double dr6, dr2;

        //cout << "neighbour size: " << atom1->neighbourList().size() << endl;
        //cout << "Neighbour size: " << m_cellList->getNeighbours()[i].size() << endl;
        // loop over all atoms in atom1's neighbour list
        //for (int j=0; j < atom1->neighbourList().size(); j++) {
        for (int j=0; j < m_cellList->getNeighbours()[i].size(); j++) {

            //Atom *atom2 = atom1->neighbourList()[j];
            Atom *atom2 = m_cellList->getNeighbours()[i][j];

            // calculate distance vector
            dr = atom1->position - atom2->position;

            // make sure we're using shortest distance component-wise (periodic boundary conditions)
            for (int dim=0; dim < 3; dim++) {
                if (dr[dim] > m_system.systemSizeHalf()[dim]) { dr[dim] -= m_system.systemSize()[dim]; }
                else if (dr[dim] < -m_system.systemSizeHalf()[dim]) { dr[dim] += m_system.systemSize()[dim]; }
            }

            // calculate force
            dr6 = 1.0 / pow(dr.length(), 6);
            dr2 = 1.0 / dr.lengthSquared();
            forceOnAtom = 24*dr6*(2*dr6 - 1)*dr2*dr;

            //cout << forceOnAtom << endl;
            // add contribution to force on atom i and j
            atom1->force += forceOnAtom;
            atom2->force -= forceOnAtom;   // Newton's third law

            // calculate pressure
            vec3 Fij = atom1->force - atom2->force;

            // dot product of Fij and dr
            pressure += Fij.dot(dr);

            // calculate potential energy
            potentialEnergy += (dr6 - 1)*dr6;
        }
    }


    m_potentialEnergy = 4*potentialEnergy;
    m_pressure = m_inverseVolume*pressure;

}
