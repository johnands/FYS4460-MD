#include "celllist.h"
#include "atom.h"
#include "system.h"
#include "math/vec3.h"
#include <iostream>

using std::cout;
using std::endl;
using std::vector;

CellList::CellList(System *system, double rCut, double neighbourCut) {
    m_system = system;
    m_rCut = rCut;
    m_rCutSquared = rCut*rCut;
    m_neighbourCut = neighbourCut;
    m_neighbourCutSquared = neighbourCut*neighbourCut;
}

void CellList::setupCells() {

    m_numberOfCellsEachDimension = m_system->systemSize().x() / m_rCut;
    if (m_numberOfCellsEachDimension < 1.0) { m_numberOfCellsEachDimension = 1.0; }

    cout << "System size: " << m_system->systemSize().x() << endl;
    cout << "Cutoff: " << m_rCut << endl;
    cout << "Number of neighbour cells: " << m_numberOfCellsEachDimension*3 << endl;

    // make room in memory
    m_cells.resize( m_numberOfCellsEachDimension,
                    vector<vector<vector<Atom*>>>(m_numberOfCellsEachDimension,
                    vector<vector<Atom*>> (m_numberOfCellsEachDimension,
                    vector<Atom*> ())) );

    m_neighbours.resize(m_system->atoms().size());
}


void CellList::updateCells() {
    // set up cell lists as a four-dimensional std vector
    vec3 systemSize = m_system->systemSize();

    for (Atom *atom : m_system->atoms()) {

        // find which cell the atom is in
        int i = atom->position.x() / systemSize.x() * m_numberOfCellsEachDimension;
        int j = atom->position.y() / systemSize.y() * m_numberOfCellsEachDimension;
        int k = atom->position.z() / systemSize.z() * m_numberOfCellsEachDimension;

        /*
        cout << "i: " << i << endl;
        cout << "j: " << j << endl;
        cout << "k: " << k << endl;
        */

        // set cell indicies
        atom->setCellIndicies(i, j, k);

        // store atom in the designated cell
        m_cells[i][j][k].push_back(atom);
    }

    // make neighbour list for each cell
    for (int ci=0; ci < m_numberOfCellsEachDimension; ci++) {
        for (int cj=0; cj < m_numberOfCellsEachDimension; cj++) {
            for (int ck=0; ck < m_numberOfCellsEachDimension; ck++) {
                for (int di=0; di <= 1; di++) {
                    for (int dj=(di==0 ? 0 : -1); dj <= 1; dj++) {
                        for (int dk=(di==0 && dj==0 ? 0 : -1); dk <= 1; dk++) {

                            // take care of the periodic boundary conditions
                            // eg. if atom in end cell, then neighbour is first cell for a given dimension
                            const int cii = (ci+di==-1 ? m_numberOfCellsEachDimension-1 : (ci+di==m_numberOfCellsEachDimension ? 0 : ci+di));
                            const int cjj = (cj+dj==-1 ? m_numberOfCellsEachDimension-1 : (cj+dj==m_numberOfCellsEachDimension ? 0 : cj+dj));
                            const int ckk = (ck+dk==-1 ? m_numberOfCellsEachDimension-1 : (ck+dk==m_numberOfCellsEachDimension ? 0 : ck+dk));

                            // choose another atom
                            // a runs through all atom objects in cell [ci, cj, ck]
                            // b runs through all atom objects in the chosen neighbour cell [cii, cjj, ckk]
                            for (int a=0; a < m_cells[ci][cj][ck].size(); a++) {
                                for (int b=(di==0 && dj==0 && dk==0 ? a+1 : 0); b < m_cells[cii][cjj][ckk].size(); b++) {

                                    Atom* atom1 = m_cells[ci][cj][ck][a];       // atom in current cell
                                    Atom* atom2 = m_cells[cii][cjj][ckk][b];    // atom in chosen neighbour cell

                                    // find distance between atom1 and atom2
                                    vec3 dr;
                                    dr = atom2->position - atom1->position;

                                    // apply periodic boundary conditions
                                    for (int dim=0; dim < 3; dim++) {
                                        if (dr[dim] > m_system->systemSizeHalf()[dim]) { dr[dim] -= m_system->systemSize()[dim]; }
                                        else if (dr[dim] < -m_system->systemSizeHalf()[dim]) { dr[dim] += m_system->systemSize()[dim]; }
                                    }

                                    // add atom2 to atom1 neighbour list if distance is shorter than cut-off distance
                                    if (dr.lengthSquared() < m_neighbourCutSquared) {
                                        m_neighbours[atom1->getIndex()].push_back(atom2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void CellList::clearCells() {
    for (int i=0; i < m_numberOfCellsEachDimension; i++) {
        for (int j=0; j < m_numberOfCellsEachDimension; j++) {
            for (int k=0; k < m_numberOfCellsEachDimension; k++) {
                /*for (int l=0; l < m_cells[i][j][k].size(); l++) {
                    m_cells[i][j][k][l]->neighbourList().clear();
                }*/
                m_cells[i][j][k].clear();
            }  
        }
    }
    /*for (int l=0; l < m_system->atoms().size(); l++) {
        m_system->atoms()[l]->neighbourList().clear();*/
    for (int l = 0; l < m_system->atoms().size(); l++) {
        m_neighbours[l].clear();
    }
}
