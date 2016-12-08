//#include "tensorflownetwork.h"
//#include "../celllist.h"
//#include "tensorflow/core/platform/env.h"
//#include <vector>

//TensorFlowNetwork::TensorFlowNetwork(System *system, std::string filename, double rCut, double neighbourCut) :
//                       Potential(system) {

//    m_filename = filename;

//    m_cellList = new CellList(system, rCut, neighbourCut);
//    m_cellList->setupCells();
//    m_cellList->updateCells();

//    m_rCut = rCut;
//    m_rCutSquared = rCut*rCut;
//    m_neighbourCut = neighbourCut;

//    double r2 = 1.0 / m_rCutSquared;
//    m_potentialCut = 4*r2*r2*r2*(r2*r2*r2 - 1);

//    m_updateLists = 0;

//    setUpNetwork();
//}


//void TensorFlowNetwork::setUpNetwork() {

//    // Initialize a tensorflow session
//    m_status = tensorflow::NewSession(tensorflow::SessionOptions(), &m_session);
//    if ( !m_status.ok() ) {
//      std::cout << m_status.ToString() << "\n";
//      std::cout << "Failed to launch session" << std::endl;
//      exit(1);
//    }

//    // Read in the protobuf graph we exported
//    tensorflow::GraphDef graph_def;
//    m_status = tensorflow::ReadBinaryProto(tensorflow::Env::Default(), m_filename,
//                                           &graph_def);
//    if ( !m_status.ok() ) {
//      std::cout << m_status.ToString() << "\n";
//      std::cout << "Failed to read proto file" << std::endl;
//      exit(2);
//    }


//    // Add the graph to the session
//    m_status = m_session->Create(graph_def);
//    if ( !m_status.ok() ) {
//      std::cout << m_status.ToString() << "\n";
//      std::cout << "Failed to add graph to session" << std::endl;
//      exit(3);
//    }
//}


//void TensorFlowNetwork::calculateForces() {

//    double potentialEnergy = 0;
//    double pressure = 0;

//    // update lists every 5th time step
//    if (m_updateLists >= 20) {
//        m_cellList->clearCells();
//        m_cellList->updateCells();
//        m_updateLists = 0;
//    }

//    m_updateLists++;

//    for (int i=0; i < m_system->atoms().size(); i++) {
//        Atom *atom1 = m_system->atoms()[i];

//        // loop over all atoms in atom1's neighbour list
//        for (int j=0; j < m_cellList->getNeighbours()[i].size(); j++) {

//            Atom *atom2 = m_cellList->getNeighbours()[i][j];

//            // calculate distance vector
//            double dr[] = {atom1->position[0] - atom2->position[0],
//                           atom1->position[1] - atom2->position[1],
//                           atom1->position[2] - atom2->position[2]};

//            // make sure we're using shortest distance component-wise (periodic boundary conditions)
//            for (int dim=0; dim < 3; dim++) {
//                if      (dr[dim] > m_system->systemSizeHalf()[dim])  { dr[dim] -= m_system->systemSize()[dim]; }
//                else if (dr[dim] < -m_system->systemSizeHalf()[dim]) { dr[dim] += m_system->systemSize()[dim]; }
//            }

//            double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

//            // calculate potential energy with ANN using distance as input if
//            // distance is shorter than cutoff
//            if (dr2 <= m_rCutSquared) {

//                // calculate distance
//                float distance = sqrt(dr2);

//                // initialize tensor to hold distance value
//                // type: float, shape: no shape, but one element
//                tensorflow::Tensor r_ij(tensorflow::DT_FLOAT, tensorflow::TensorShape({1,1}));

//                // represent distance as a tensorflow scalar tensor
//                r_ij.scalar<float>()() = distance;

//                // input to Run method needs to be a std::vector of std::pairs
//                std::vector<std::pair<string, tensorflow::Tensor>> inputs = {
//                  { "input/x-input:0", r_ij },
//                };

//                // the session will initialize the outputs
//                std::vector<tensorflow::Tensor> outputs;

//                clock_t begin = clock();
//                // run the session, evaluating the network with given inputs
//                m_status = m_session->Run(inputs, {"outputLayer/activation"}, {}, &outputs);
//                clock_t end = clock();
//                double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//                std::cout << "Did som energies yeah " << elapsed_secs << std::endl;
//                if (!m_status.ok()) {
//                  std::cout << m_status.ToString() << "\n";
//                  std::cout << "Failed to evalute network" << std::endl;
//                  exit(4);
//                }

//                // grab the first output (we only evaluated one graph node: "prediction")
//                // and convert the node to a scalar representation.
//                auto energy = outputs[0].scalar<float>();
//                //potentialEnergy += energy;
//                // forceOnAtom = dr*0;

//                // Print the results
//                //std::cout << outputs[0].DebugString() << "\n";
//                //std::cout << energy << "\n";

//                // Free any resources used by the session
//                //m_session->Close();

//            }
//            else {
//                // if not, zero force and energy
//                // forceOnAtom = dr*0;
//            }

//            // add contribution to force on atom i and j
////            atom1->force += forceOnAtom;
////            atom2->force -= forceOnAtom;   // Newton's third law

//            if (m_system->getUseExternalForce()) {
//                atom1->force[0] += 1.0;
//            }

//            // dot product of Fij and dr
//            // pressure += forceOnAtom.dot(dr);
//        }
//    }

//    m_potentialEnergy = potentialEnergy;
//    m_pressure = m_inverseVolume*pressure;
//}
