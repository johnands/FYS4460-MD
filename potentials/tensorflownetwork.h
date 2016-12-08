//#ifndef TENSORFLOWNETWORK_H
//#define TENSORFLOWNETWORK_H
//#include "potential.h"
//#include "tensorflow/core/public/session.h"
//#include <string>

//class CellList;

//class TensorFlowNetwork : public Potential {

//public:
//    TensorFlowNetwork(System *system, std::string filename, double rCut, double neighbourCut);
//    void calculateForces();
//    void setUpNetwork();

//private:
//    CellList *m_cellList = nullptr;
//    int    m_updateLists = 0;
//    double m_rCut = 0.0;
//    double m_rCutSquared = 0.0;
//    double m_neighbourCut = 0.0;
//    double m_potentialCut = 0.0;

//    std::string m_filename;
//    tensorflow::Session *m_session = nullptr;
//    tensorflow::Status m_status;


//};

//#endif // TENSORFLOWNETWORK_H
