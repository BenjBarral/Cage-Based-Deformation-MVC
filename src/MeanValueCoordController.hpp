//
//  MeanValueCoordController.hpp
//  3DA_project_CageBasedDef_bin
//
//  Created by Benjamin Barral on 11/05/2019.
//

#ifndef MeanValueCoordController_hpp
#define MeanValueCoordController_hpp

#include <igl/readPLY.h>
#include <cmath>

using namespace Eigen;
using namespace std;

class MeanValueCoordController{
private:
    MatrixXd V_mesh, V_cage,V_cage_deformed;
    MatrixXi F_mesh, F_cage;
    MatrixXd mV_weights;
    int num_vertices_mesh,num_vertices_cage;
    int num_faces_mesh,num_faces_cage;
    double epsilon_;
    
public:
    MeanValueCoordController();
    MeanValueCoordController(const MatrixXd& vertices_mesh, const MatrixXd& vertices_cage,
                             const MatrixXi& faces_mesh, const MatrixXi& faces_cage,
                             const double& bb_size);
    void ComputeMVWeights();
    void SetDeformedCage(const MatrixXd& vertices_cage_deformed);
    MatrixXd MVInterpolate();
    
};

#endif /* MeanValueCoordController_hpp */
