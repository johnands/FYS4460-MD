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
    double m_rCut = 0.0;
    double m_rCutSquared = 0.0;
    double m_neighbourCut = 0.0;
    double m_potentialCut = 0.0;
    int    m_updateLists = 0;
    CellList *m_cellList = nullptr;

public:
    LennardJonesCellList(System *system, double sigma, double epsilon, double rCut, double neighbourCut);
    virtual void calculateForces();

    // setters and getters
    CellList *celllist() { return m_cellList; }
    void setCellList(CellList *celllist) { m_cellList = celllist; }
};

#endif
