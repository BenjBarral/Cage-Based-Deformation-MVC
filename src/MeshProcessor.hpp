//
//  MeshProcessor.hpp
//  3DA_project_CageBasedDef_bin
//
//  Created by Benjamin Barral on 12/05/2019.
//

#ifndef MeshProcessor_hpp
#define MeshProcessor_hpp

#include <stdio.h>
#include <igl/readPLY.h>

// A class that computes the barycenter and extreme coordinates of a mesh

using namespace Eigen;
using namespace std;

struct ExtremeVertex{
    int vertex_index_;
    double coord_;
    ExtremeVertex(){
    }
};

class MeshProcessor{
private:
    MatrixXd V_;
    MatrixXi F_;
    int num_vertices_, num_faces_;
    Vector3d barycenter_;
    ExtremeVertex* max_vertices_, *min_vertices_;
    
    
public:
    MeshProcessor();
    MeshProcessor(const MatrixXd& vertices, const MatrixXi& faces);
    Vector3d GetBarycenter();
    void ComputeExtremeVertices();
    int GetMaximumVertexIndex(const int& dim);
    int GetMinimumVertexIndex(const int& dim);
    double GetBoundingBoxSize();
    
};
#endif /* MeshProcessor_hpp */
