#ifndef LENNARDJONESCELLLIST_H
#define LENNARDJONESCELLLIST_H
#include "potential.h"

class CellList;

class LennardJonesCellList : public Potential
{
private:
    double m_sigma   = 1.0;
    double m_epsilon = 1.0;
    double m_sigma6  = 1.0;
    double m_cutOffDistance = 0.0;
    CellList *m_cellList = nullptr;

public:
    LennardJonesCellList(System &system, double sigma, double epsilon, double cutOffDistance);
    virtual void calculateForces();

    // setters and getters
    CellList *celllist() { return m_cellList; }
    void setCellList(CellList *celllist) { m_cellList = celllist; }
};

#endif
