#ifndef CELLLIST_H
#define CELLLIST_H
#include <vector>
#include "system.h"

using std::vector;

class CellList {

public:
    CellList(class System &system, double cutOffDistance);
    void setupCells();
    void updateCells();
    void clearCells();
    vector<class Atom *>& getCell(int i, int j, int k) { return m_cells[i][j][k]; }

private:
    vector<vector<vector<vector<class Atom *>>>> m_cells;
    int m_numberOfCellsEachDimension;
    double m_cutOffDistance;
    double m_cutOffDistance2;
    class System m_system;
};

#endif // CELLLIST_H
