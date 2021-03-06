TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXX += -g

release {
    DEFINES += ARMA_NO_DEBUG
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}

SOURCES += main.cpp \
    atom.cpp \
    system.cpp \
    integrators/integrator.cpp \
    integrators/velocityverlet.cpp \
    math/vec3.cpp \
    math/random.cpp \
    io.cpp \
    potentials/potential.cpp \
    potentials/lennardjones.cpp \
    statisticssampler.cpp \
    integrators/eulercromer.cpp \
    unitconverter.cpp \
    celllist.cpp \
    potentials/lennardjonescelllist.cpp \
    thermostats/thermostat.cpp \
    thermostats/berendsen.cpp \
    thermostats/andersen.cpp \
    porosities/porosities.cpp \
    porosities/centeredcylinder.cpp \
    porosities/spheres.cpp \
    potentials/neuralnetwork.cpp \
    math/activationfunctions.cpp \
    potentials/tensorflownetwork.cpp \
    examples.cpp \
    potentials/manyneighbournn.cpp

HEADERS += \
    atom.h \
    system.h \
    integrators/integrator.h \
    integrators/velocityverlet.h \
    math/vec3.h \
    math/random.h \
    io.h \
    potentials/potential.h \
    potentials/lennardjones.h \
    statisticssampler.h \
    integrators/eulercromer.h \
    unitconverter.h \
    celllist.h \
    potentials/lennardjonescelllist.h \
    thermostats/thermostat.h \
    thermostats/berendsen.h \
    thermostats/andersen.h \
    porosities/porosities.h \
    porosities/centeredcylinder.h \
    porosities/spheres.h \
    potentials/neuralnetwork.h \
    math/activationfunctions.h \
    potentials/tensorflownetwork.h \
    examples.h \
    potentials/manyneighbournn.h

#INCLUDEPATH += /home/johnands/QTensorFlow/include
#LIBS += -L/home/johnands/QTensorFlow/lib64 -ltensorflow

LIBS += -larmadillo -lblas -llapack

