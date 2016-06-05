#include "io.h"
#include "system.h"
#include "atom.h"
#include "unitconverter.h"
#include <cstdlib>

using std::endl; using std::cout;

IO::IO()
{

}

IO::~IO() {
    close();
}

void IO::open(const char *filename) {
    if(file.is_open()) {
        std::cout << "<IO.cpp> Error, tried to open file " << filename << ", but some file is already open." << endl;
        exit(1);
    }

    file.open(filename);
}

void IO::close() {
    if(file.is_open()) {
        file.close();
    }
}

// This saves the current state to a file following the xyz-standard (see http://en.wikipedia.org/wiki/XYZ_file_format )
void IO::saveState(System *system)
{
    if(file.is_open()) {
        if (system->steps() == 501) {
        file << system->atoms().size() << endl;
        file << "x y z vx vy vz" << endl;
        for(Atom *atom : system->atoms()) {
            file << atom->getName() << " "
                 << UnitConverter::lengthToAngstroms(atom->position.x()) << " " << UnitConverter::lengthToAngstroms(atom->position.y()) << " " << UnitConverter::lengthToAngstroms(atom->position.z()) << " "
                 << atom->velocity.x() << " " << atom->velocity.y() << " " << atom->velocity.z() <<  endl;
        }}
    }
}
