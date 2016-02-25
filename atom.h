#ifndef ATOM_H
#define ATOM_H
#include "math/vec3.h"
#include <vector>

class Atom
{
private:
    float m_mass;
    std::vector<Atom *> m_neighbourList;
    vec3 m_cellIndicies;

public:
    vec3 position;
    vec3 velocity;
    vec3 force;

    Atom(double mass);
    void resetForce();
    void resetVelocityMaxwellian(double temperature);

    // getters
    double mass() { return m_mass; }
    std::vector<Atom *> neighbourList() { return m_neighbourList; }
    vec3 cellIndicies() { return m_cellIndicies; }

    // setters
    void setMass(double mass) { m_mass = mass; }
    void addToNeighbourList(Atom * atom) { m_neighbourList.push_back(atom); }
    void setCellIndicies(int i, int j, int k);
};
#endif
