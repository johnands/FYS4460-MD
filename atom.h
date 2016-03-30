#ifndef ATOM_H
#define ATOM_H
#include "math/vec3.h"
#include <vector>

class Atom
{
private:
    float m_mass;
    int m_index;
    std::vector<Atom *> m_neighbourList;
    vec3 m_cellIndicies;
    vec3 m_initialPosition;
    bool m_movingAtom = true;
    const char *m_name = "Ar";

public:
    vec3 position;
    vec3 velocity;
    vec3 force;

    Atom(double mass);
    void resetForce();
    void resetVelocityMaxwellian(double temperature);
    void resetVelocityUniform(double maxMinVelocity);
    void storeInitialPosition();

    // getters
    double mass() { return m_mass; }
    int getIndex() { return m_index; }
    std::vector<Atom *> neighbourList() { return m_neighbourList; }
    vec3 cellIndicies() { return m_cellIndicies; }
    vec3 initialPosition() { return m_initialPosition; }
    bool movingAtom() { return m_movingAtom; }
    const char *getName() { return m_name; }

    // setters
    void setMass(double mass) { m_mass = mass; }
    void addToNeighbourList(Atom * atom) { m_neighbourList.push_back(atom); }
    void setCellIndicies(int i, int j, int k);
    void setIndex(int index);
    void adjustIndex(int change);
    void setMovingAtom(bool moving);
    void setName(const char *name);
};
#endif
