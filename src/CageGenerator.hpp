//
//  CageGenerator.hpp
//  3DA_project_CageBasedDef_bin
//
//  Created by Benjamin Barral on 13/05/2019.
//

// Automated coarse cage generator based on https://www.cse.wustl.edu/~taoju/research/meanvalue.pdf

#ifndef CageGenerator_hpp
#define CageGenerator_hpp

#include <stdio.h>
#include <igl/readPLY.h>
#include <Eigen/QR>
#include <Eigen/Sparse>
#include <igl/cotmatrix.h>
#include "igl/ray_mesh_intersect.h"
//#include "igl/copyleft/cgal/intersect_other.h"

using namespace Eigen;
using namespace std;

class CageGenerator{
private:
    double sparseness_;
    MatrixXd V_mesh_,V_pca_basis_, V_cage_, V_cage_smooth_;
    MatrixXi F_mesh_, F_cage_;
    int num_vertices_mesh_, num_faces_mesh_;
    int num_vertices_cage_, num_faces_cage_;
    Vector3d barycenter_;
    MatrixXd pca_basis_matrix_;
    Vector3d bbBox_ptMin, bbBox_ptMax;
    int* n_vox;
    double* res_vox,* epsilon_tests;
    vector < vector < vector<Vector3d> > > voxel_center_positions_, voxel_face_vertices_;
    MatrixXd grid_face_vertices_;
    MatrixXd voxel_face_normals_;
    vector < vector < vector<int> > > voxel_partition_info_;
    MatrixXd outer_mesh_vertices_;
    MatrixXi outer_mesh_faces_;
    SparseMatrix<double> laplace_beltrami_matrix_;
    int num_smoothing_iterations_;
    double lambda_smooth_;
    
public:
    CageGenerator(const MatrixXd& vertices, const MatrixXi& faces, const int& smooth_it, const double& lambda,  const double& sparseness);
    void ComputePCA();
    void ComputePrincipalDirectionsBoundingBox(MatrixXd& bb_vertices, MatrixXi& bb_faces);
    void ComputeCubeMesh(MatrixXd& bb_vertices, MatrixXi& bb_faces, const Vector3d& min_point, const Vector3d& max_point);
    bool IsPtInsideVoxel(const Vector3d& point, const Vector3d& vox_origin);
    bool IsPtOnVoxelSurface(const Vector3d& point, const Vector3d& vox_origin);
    bool IsPtInsideBoudingBox(const Vector3d& point);
    bool VoxelTriangleIntersection(const MatrixXd& triangle, const Vector3d& vox_origin);
    void InitializeVoxelPositions();
    void InitializeVoxelNormals();
    int TupleToIndex(const int& i, const int& j, const int& k);
    bool VoxelMeshIntersection(const Vector3d& vox_orig);
    void FindFeatureVoxels(vector<MatrixXd>& feat_voxels_vertices, vector<MatrixXi>& feat_voxels_faces);
    void ComputeSeparatingFaces(const Vector3i& curr_vox_tuple, const int& q,
                                Vector3i& face_1, Vector3i& face_2, const bool& reverse);
    void ExtractOuterSurface();
    void SimplifyAdjacentFaces();
    void FillingAlgorithm(const Vector3i& vox_ind, vector<Vector3i>& outer_faces);
    void ComputeCage(MatrixXd& bb_vertices, MatrixXi& bb_faces,
                     vector<MatrixXd>& feat_voxels_vertices, vector<MatrixXi>& feat_voxels_faces);
    MatrixXd GetGridVertices();
    void SmoothCage();
    MatrixXd GetCage();
    void SetSmoothingParameters(const int& n_i, const float& lambda);
    MatrixXd GetSmoothCage();
    MatrixXi GetCageFaces();
};

#endif /* CageGenerator_hpp */
