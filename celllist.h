#ifndef CELLLIST_H
#define CELLLIST_H
#include <vector>
#include "system.h"

using std::vector;

class CellList {

public:
    CellList(class System *system, double rCut, double neighbourCut);
    void setupCells();
    void updateCells();
    void clearCells();
    vector<class Atom *>& getCell(int i, int j, int k) { return m_cells[i][j][k]; }
    vector<vector<class Atom *>>& getNeighbours() { return m_neighbours; }

private:
    vector<vector<vector<vector<class Atom *>>>> m_cells;
    vector<vector<class Atom*>> m_neighbours;
    int m_numberOfCellsEachDimension;
    double m_rCut;
    double m_rCutSquared;
    double m_neighbourCut;
    double m_neighbourCutSquared;
    class System *m_system = nullptr;
};

#endif // CELLLIST_H
