#include "celllist.h"
#include "atom.h"
#include "system.h"
#include "math/vec3.h"
#include <iostream>

using std::cout;
using std::endl;

CellList::CellList(System *system, double cutOffDistance) {
    m_system = system;
    m_cutOffDistance = cutOffDistance;
}

void CellList::setup() {
    // set up cell lists as a four-dimensional std vector
    vec3 systemSize = m_system->systemSize();

    m_numberOfCellsEachDimension = systemSize.x() / m_cutOffDistance;

    for (Atom *atom : m_system->atoms()) {

        // find which cell the atom is in
        int i = atom->position.x() / systemSize.x() * m_numberOfCellsEachDimension;
        int j = atom->position.y() / systemSize.y() * m_numberOfCellsEachDimension;
        int k = atom->position.z() / systemSize.z() * m_numberOfCellsEachDimension;

        cout << i << endl << j << endl << k << endl;

        // set cell indicies
        atom->setCellIndicies(i, j, k);

        // store atom in the designated cell
        m_cells[i][j][k].push_back(atom);
    }
}

CellList::clearCells() {
    for (int i=0; i < m_numberOfCellsEachDimension; i++) {
        for (int j=0; j < m_numberOfCellsEachDimension; j++) {
            for (int k=0; k < m_numberOfCellsEachDimension; k++) {
                m_cells[i][j][k].clear();
            }
        }
    }
}
