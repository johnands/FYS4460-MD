#ifndef CELLLIST_H
#define CELLLIST_H
#include <vector>

using std::vector;
class System;

class CellList {

public:
    CellList(System *system, double cutOffDistance);
    void setup();
    void clearCells();
    vector<class Atom *>& getCell(int i, int j, int k) { return m_cells[i][j][k]; }

private:
    vector<vector<vector<vector<class Atom *>>>> m_cells;
    int m_numberOfCellsEachDimension;
    double m_cutOffDistance;
    System *m_system;
};

#endif // CELLLIST_H
