//
//  CageGenerator.cpp
//  3DA_project_CageBasedDef_bin
//
//  Created by Benjamin Barral on 13/05/2019.
//

// Automated coarse cage generator based on https://www.cse.wustl.edu/~taoju/research/meanvalue.pdf

#include "CageGenerator.hpp"

#define NUM_FACES_CUBE 6
#define NUM_FACES_CUBE_TRI_MESH 12
#define NUM_VERTICES_CUBE 8

static const double EPSILON_TESTS_RATIO = 0.001;

static bool implicit_smoothing = true;
static const double EPSILON_BOUNDING_BOX_ = 0.15;  // make the bounding box a bit larger

CageGenerator::CageGenerator(const MatrixXd& vertices, const MatrixXi& faces, const int& smooth_it, const double& lambda, const double& sparseness) {
    sparseness_ = sparseness;
    V_mesh_ = vertices;
    F_mesh_ = faces;
    num_vertices_mesh_ = vertices.rows();
    num_faces_mesh_ = faces.rows();
    num_smoothing_iterations_ = smooth_it;
    
    barycenter_ = Vector3d::Zero(3);
    for (int i = 0; i < num_vertices_mesh_; i++){
        barycenter_ += vertices.row(i);
    }
    barycenter_ /= num_vertices_mesh_;
    InitializeVoxelNormals();
    lambda_smooth_ = lambda;
}

void CageGenerator::InitializeVoxelNormals() {
    voxel_face_normals_ = MatrixXd::Zero(NUM_FACES_CUBE,3);
    voxel_face_normals_.row(0) = Vector3d::UnitX();
    voxel_face_normals_.row(1) = -Vector3d::UnitX();
    voxel_face_normals_.row(2) = Vector3d::UnitY();
    voxel_face_normals_.row(3) = -Vector3d::UnitY();
    voxel_face_normals_.row(4) = Vector3d::UnitZ();
    voxel_face_normals_.row(5) = -Vector3d::UnitZ();
}

void CageGenerator::InitializeVoxelPositions() {
    voxel_center_positions_.clear();
    voxel_face_vertices_.clear();
    Vector3d deltaVox(res_vox[0]/2., res_vox[1]/2., res_vox[2]/2.);
    double x,y,z;
    Vector3d vox_orig,vox_vertex;
    for (int i = 0; i < n_vox[0]; i++){
        voxel_center_positions_.push_back(vector<vector<Vector3d>>());
        x = bbBox_ptMin(0) + ( ((double)i + 0.5 )/ n_vox[0] ) * (bbBox_ptMax(0) - bbBox_ptMin(0));
        for (int j = 0; j < n_vox[1]; j++){
            voxel_center_positions_[i].push_back(vector<Vector3d>());
            y = bbBox_ptMin(1) + ( ((double)j + 0.5 )/ n_vox[1] ) * (bbBox_ptMax(1) - bbBox_ptMin(1));
            for (int k = 0; k < n_vox[2]; k++){
                z = bbBox_ptMin(2) + ( ((double)k + 0.5 )/ n_vox[2] ) * (bbBox_ptMax(2) - bbBox_ptMin(2));
                vox_orig = Vector3d(x,y,z);
                voxel_center_positions_[i][j].push_back(vox_orig);
            }
        }
    }
    int ind;
    int num_vert_grid = (n_vox[0]+1) * (n_vox[1]+1) * (n_vox[2]+1);
    grid_face_vertices_ = MatrixXd::Zero(num_vert_grid, 3);
    for (int i = 0; i <= n_vox[0]; i++){
        voxel_face_vertices_.push_back(vector<vector<Vector3d>>());
        x = bbBox_ptMin(0) + ( ((double)i)/ n_vox[0] ) * (bbBox_ptMax(0) - bbBox_ptMin(0));
        for (int j = 0; j <= n_vox[1]; j++){
            voxel_face_vertices_[i].push_back(vector<Vector3d>());
            y = bbBox_ptMin(1) + ( ((double)j )/ n_vox[1] ) * (bbBox_ptMax(1) - bbBox_ptMin(1));
            for (int k = 0; k <= n_vox[2]; k++){
                z = bbBox_ptMin(2) + ( ((double)k)/ n_vox[2] ) * (bbBox_ptMax(2) - bbBox_ptMin(2));
                vox_vertex = Vector3d(x,y,z);
                voxel_face_vertices_[i][j].push_back(vox_vertex);
                ind = TupleToIndex(i, j, k);
                grid_face_vertices_.row(ind) = vox_vertex;
            }
        }
    }
    num_vertices_cage_ = num_vert_grid;
}

void CageGenerator::ComputeCubeMesh(MatrixXd &bb_vertices, MatrixXi &bb_faces, const Vector3d &min_point, const Vector3d &max_point) {
    int n_f = NUM_FACES_CUBE_TRI_MESH;
    bb_vertices = MatrixXd::Zero(NUM_VERTICES_CUBE, 3);
    bb_faces = MatrixXi::Zero(n_f, 3);
    // Faces
    bb_faces.row(0) = Vector3i(0,1,2);
    bb_faces.row(1) = Vector3i(0,2,3);
    bb_faces.row(2) = Vector3i(1,5,6);
    bb_faces.row(3) = Vector3i(1,6,2);
    bb_faces.row(4) = Vector3i(5,4,7);
    bb_faces.row(5) = Vector3i(5,7,6);
    bb_faces.row(6) = Vector3i(4,0,3);
    bb_faces.row(7) = Vector3i(4,3,7);
    bb_faces.row(8) = Vector3i(4,5,1);
    bb_faces.row(9) = Vector3i(4,1,0);
    bb_faces.row(10) = Vector3i(3,2,6);
    bb_faces.row(11) = Vector3i(3,6,7);
    
    // Vertices
    double xMin = min_point[0], yMin = min_point[1], zMin = min_point[2];
    double xMax = max_point[0], yMax = max_point[1], zMax = max_point[2];
    bb_vertices.row(0) = Vector3d(xMin, yMin, zMin);
    bb_vertices.row(1) = Vector3d(xMax, yMin, zMin);
    bb_vertices.row(2) = Vector3d(xMax, yMin, zMax);
    bb_vertices.row(3) = Vector3d(xMin, yMin, zMax);
    bb_vertices.row(4) = Vector3d(xMin, yMax, zMin);
    bb_vertices.row(5) = Vector3d(xMax, yMax, zMin);
    bb_vertices.row(6) = Vector3d(xMax, yMax, zMax);
    bb_vertices.row(7) = Vector3d(xMin, yMax, zMax);
}

void CageGenerator::ComputePCA() { 
    MatrixXd V_bar = V_mesh_;
    pca_basis_matrix_ = MatrixXd::Zero(3, 3);
    for (int i = 0; i < num_vertices_mesh_; i++){
        V_bar.row(i) -= barycenter_;
    }
    MatrixXd cov_mat = V_bar.transpose() * V_bar / num_vertices_mesh_;
    EigenSolver<MatrixXd> eigen_solver(cov_mat);
    for (int i = 0; i < 3; i++){
        Vector3d eigVec = eigen_solver.eigenvectors().col(i).real();
        pca_basis_matrix_.col(i) = eigVec / eigVec.norm();
    }
    V_pca_basis_ = (pca_basis_matrix_.transpose() * V_bar.transpose()).transpose();
}

void CageGenerator::ComputePrincipalDirectionsBoundingBox(MatrixXd &bb_vertices, MatrixXi &bb_faces) {
    ComputePCA();
    // Change basis
    res_vox = new double[3];
    n_vox = new int[3];
    double minCoords[3] = {__DBL_MAX__, __DBL_MAX__, __DBL_MAX__};
    double maxCoords[3] = {__DBL_MIN__, __DBL_MIN__, __DBL_MIN__};
    Vector3d point;
    for (int i = 0; i < num_vertices_mesh_; i++){
        point = V_pca_basis_.row(i);
        for (int j = 0; j < 3; j++){
            if (point(j) > maxCoords[j]){
                maxCoords[j] = point(j);
            }
            if (point(j) < minCoords[j]){
                minCoords[j] = point(j);
            }
        }
    }
    bbBox_ptMin = Vector3d::Zero(3);
    bbBox_ptMax = Vector3d::Zero(3);
    epsilon_tests = new double[3];
    Vector3d deltaBB(maxCoords[0] - minCoords[0], maxCoords[1] - minCoords[1], maxCoords[2] - minCoords[2]);
    int j_max;
    double dim_max = __DBL_MIN__;
    for (int j = 0; j < 3; j++){
        bbBox_ptMin(j) = minCoords[j] - EPSILON_BOUNDING_BOX_ * deltaBB(j);
        bbBox_ptMax(j) = maxCoords[j] + EPSILON_BOUNDING_BOX_ * deltaBB(j);
        double dim_size = bbBox_ptMax(j) - bbBox_ptMin(j);
        if (dim_size > dim_max){
            dim_max = dim_size;
            j_max = j;
        }
    }
    
    n_vox[j_max] = (int)(sqrt(sparseness_ * double(num_vertices_mesh_) / 6.)); // Chuhua Xian. https://www.cse.wustl.edu/~taoju/research/meanvalue.pdf
    cout << "n_vox = " << n_vox[j_max] << endl;
    res_vox[j_max] = (bbBox_ptMax(j_max) - bbBox_ptMin(j_max)) / n_vox[j_max];
    for (int j = 0; j < 3; j++){
        if (j != j_max){
            int n = (int)((bbBox_ptMax(j) - bbBox_ptMin(j)) / res_vox[j_max]) + 2;
            n_vox[j] = n;
            cout << "n_vox = " << n_vox[j] << endl;
            res_vox[j] = (bbBox_ptMax(j) - bbBox_ptMin(j)) / n_vox[j];
            epsilon_tests[j] = EPSILON_TESTS_RATIO * res_vox[j];
        }
    }
    ComputeCubeMesh(bb_vertices, bb_faces, bbBox_ptMin, bbBox_ptMax);
    // Change basis
    bb_vertices = (pca_basis_matrix_ * bb_vertices.transpose()).transpose();
    for (int i = 0; i < bb_vertices.rows(); i++){
        bb_vertices.row(i) = barycenter_ +  Vector3d(bb_vertices.row(i));
    }
}

bool CageGenerator::VoxelTriangleIntersection(const MatrixXd& triangle, const Vector3d &vox_origin) {
    // If one of the three vertices is inside the voxel : true
    for (int i = 0; i < 3; i++){
        Vector3d p = triangle.row(i);
        if (IsPtInsideVoxel(p, vox_origin)) return true;
    }
    // Check intersection against the 6 voxel faces
    MatrixXd current_clipped_triangle = triangle;
    double d, sdJ, sdJPlus, tInt;
    Vector3d pJ, pJPlus, n, newPoint;
    for (int i = 0; i < NUM_FACES_CUBE; i++){
        bool allPositivePlane = true;
        n = voxel_face_normals_.row(i);
        d = - (vox_origin.dot(n) + res_vox[i/2] / 2.);
        // Test for each edge and clip
        for (int j = 0; j < 3; j++){
            pJ = current_clipped_triangle.row(j);
            pJPlus = current_clipped_triangle.row((j+1)%3);
            sdJ = pJ.dot(n) + d;
            sdJPlus = pJPlus.dot(n) + d;
            if (sdJ > 0 && sdJPlus < 0){ // enter the cube
                allPositivePlane = false;
                tInt = - (d + pJ.dot(n)) / ((pJPlus - pJ).dot(n));
                newPoint = pJ + tInt * (pJPlus - pJ);
                if (IsPtOnVoxelSurface(newPoint, vox_origin)) return true;
                current_clipped_triangle.row(j) = newPoint;
            }
            else if (sdJ < 0 && sdJPlus > 0){
                allPositivePlane = false;
                tInt = - (d + pJ.dot(n)) / ((pJPlus - pJ).dot(n));
                newPoint = pJ + tInt * (pJPlus - pJ);
                if (IsPtOnVoxelSurface(newPoint, vox_origin)) return true;
                current_clipped_triangle.row((j+1)%3) = newPoint;
            }
        }
        if (allPositivePlane) return false;
    }
    cout << "problem code" << endl;
    return false;
}

bool CageGenerator::IsPtInsideVoxel(const Vector3d &point, const Vector3d &vox_origin) {
    for (int j = 0; j < 3; j++){
        if (abs(vox_origin(j) - point(j)) > epsilon_tests[j] + res_vox[j] / 2.){
            return false;
        }
    }
    return true;
}

bool CageGenerator::IsPtOnVoxelSurface(const Vector3d& point, const Vector3d& vox_origin) {
    int count = 0;
    for (int j = 0; j < 3; j++){
        if (abs(vox_origin(j) - point(j)) > 5. * epsilon_tests[j] + res_vox[j] / 2.){
            return false;
        }
        else if (abs(vox_origin(j) - point(j)) > epsilon_tests[j] + res_vox[j] / 2.){
            count++;
        }
    }
    if (count > 1) {
        return false;
    }
    return true;
}

bool CageGenerator::VoxelMeshIntersection(const Vector3d& vox_orig) {
    Vector3d vox_min = Vector3d::Zero(3), vox_max = Vector3d::Zero(3);
    MatrixXd vox_V;
    MatrixXi vox_F;
    for (int i = 0; i < 3; i++){
        vox_min(i) = vox_orig(i) - res_vox[i] / 2.;
        vox_max(i) = vox_orig(i) + res_vox[i] / 2.;
    }
    ComputeCubeMesh(vox_V, vox_F, vox_min, vox_max);
    MatrixXi IF;
    //return intersect_other(vox_V, vox_F, V_mesh_,  F_mesh_, true, IF);
    return true;
}

void CageGenerator::FindFeatureVoxels(vector<MatrixXd>& feat_voxels_vertices, vector<MatrixXi>& feat_voxels_faces) { 
    voxel_partition_info_.clear();
    feat_voxels_vertices.clear();
    feat_voxels_faces.clear();
    int voxel_info;
    Vector3d vox_orig;
    MatrixXd triangle = MatrixXd::Zero(3,3), voxel_vertices;
    MatrixXi voxel_faces, vox_fac_tem;
    
    Vector3d deltaVox(res_vox[0]/2., res_vox[1]/2., res_vox[2]/2.);
    for (int i = 0; i < n_vox[0]; i++){
        voxel_partition_info_.push_back(vector<vector<int>>());
        for (int j = 0; j < n_vox[1]; j++){
            voxel_partition_info_[i].push_back(vector<int>());
            for (int k = 0; k < n_vox[2]; k++){
                voxel_info = -1;
                vox_orig = voxel_center_positions_[i][j][k];
                for (int f = 0; f < num_faces_mesh_; f ++){
                    for (int l = 0; l <3; l++){
                        triangle.row(l) = V_pca_basis_.row(F_mesh_(f, l));
                    }
                    if (VoxelTriangleIntersection(triangle, vox_orig)){
                        // Voxel is feature voxel
                        voxel_info = 0;
                        ComputeCubeMesh(voxel_vertices, voxel_faces, vox_orig - deltaVox, vox_orig + deltaVox);
                        // Change basis
                        voxel_vertices = (pca_basis_matrix_ * voxel_vertices.transpose()).transpose();
                        for (int i = 0; i < voxel_vertices.rows(); i++){
                            voxel_vertices.row(i) =  barycenter_ + Vector3d(voxel_vertices.row(i));
                        }
                        vox_fac_tem = voxel_faces;
                        voxel_faces.col(0) = vox_fac_tem.col(2);
                        voxel_faces.col(2) = vox_fac_tem.col(0);
                        feat_voxels_vertices.push_back(voxel_vertices);
                        feat_voxels_faces.push_back(voxel_faces);
                        break;
                    }
                }
                voxel_partition_info_[i][j].push_back(voxel_info);
            }
        }
    }
    int num_feature_voxels = feat_voxels_vertices.size();
    cout << "num_feature_voxels = " << num_feature_voxels << endl;
}

void CageGenerator::ComputeSeparatingFaces(const Vector3i& curr_vox_tuple, const int& q,
                                           Vector3i& face_1, Vector3i& face_2, const bool& reverse) {
    int i = curr_vox_tuple(0), j = curr_vox_tuple(1), k = curr_vox_tuple(2);
    int i0, i1, i2, i3;
    int j0, j1, j2, j3;
    int k0, k1, k2, k3;
    // Separating on X-axis
    if(q == 0){
        i0 = i1 = i2 = i3 = i+1;
        k0 = k1 = k;
        k2 = k3 = k + 1;
        j1 = j3 = j;
        j0 = j2 = j + 1;
    }
    else if (q==1){
        i0 = i1 = i2 = i3 = i;
        k0 = k1 = k;
        k2 = k3 = k + 1;
        j1 = j3 = j + 1;
        j0 = j2 = j;
    }
    
    // Separating on Y-axis
    if(q == 2){
        j0 = j1 = j2 = j3 = j+1;
        i0 = i2 = i;
        i1 = i3 = i+1;
        k0 = k1 = k;
        k2 = k3 = k+1;
    }
    else if (q==3){
        j0 = j1 = j2 = j3 = j;
        i0 = i2 = i+1;
        i1 = i3 = i;
        k0 = k1 = k;
        k2 = k3 = k+1;
    }
    
    // Separating on Z-axis
    if(q == 4){
        k0 = k1 = k2 = k3 = k+1;
        i0 = i2 = i;
        i1 = i3 = i+1;
        j2 = j3 = j;
        j0 = j1 = j+1;
    }
    else if (q==5){
        k0 = k1 = k2 = k3 = k;
        i0 = i2 = i+1;
        i1 = i3 = i;
        j2 = j3 = j;
        j0 = j1 = j+1;
    }
    
    face_1 = Vector3i::Zero(3), face_2 = Vector3i::Zero(3);
    int ind1 = 1, ind2 = 2;
    if (reverse){
        ind1 = 2;
        ind2 = 1;
    }
    face_1(0) = TupleToIndex(i0,j0,k0);
    face_1(ind1) = TupleToIndex(i1,j1,k1);
    face_1(ind2) = TupleToIndex(i2,j2,k2);
    face_2(0) = TupleToIndex(i1,j1,k1);
    face_2(ind1) = TupleToIndex(i3,j3,k3);
    face_2(ind2) = TupleToIndex(i2,j2,k2);
}

void CageGenerator::ExtractOuterSurface() {
    int i,j,q;
    Vector3i seed_index_tuple;
    vector<Vector3i> outer_faces;
    outer_faces.clear();
    // Expand from seeds
    for (int k = 0; k < 2; k++){
        i = k * (n_vox[0] - 1);
        for (int l = 0; l < 2; l++){
            j = l * (n_vox[1] - 1);
            for (int r = 0; r < 2; r++){
                q = r * (n_vox[2] - 1);
                seed_index_tuple = Vector3i(i,j,q);
                FillingAlgorithm(seed_index_tuple, outer_faces);
            }
        }
    }
    // Add neighboring faces from feature-to-BB
    int voxel_info;
    Vector3i curr_vox_tuple, neighbor_vox_tuple, face_1, face_2;
    bool is_neighbor_inside_bb;
    for (int i = 0; i < n_vox[0]; i++){
        for (int j = 0; j < n_vox[1]; j++){
            for (int k = 0; k < n_vox[2]; k++){
                voxel_info = voxel_partition_info_[i][j][k];
                curr_vox_tuple = Vector3i(i,j,k);
                if (voxel_info == 0){
                    for (int q = 0; q < NUM_FACES_CUBE; q++){
                        neighbor_vox_tuple = curr_vox_tuple + Vector3i(voxel_face_normals_.row(q).cast<int>());
                        is_neighbor_inside_bb = true;
                        for (int l = 0; l < 3; l++){
                            if (neighbor_vox_tuple(l) < 0 || neighbor_vox_tuple(l) > n_vox[l] - 1){
                                is_neighbor_inside_bb = false;
                                break;
                            }
                        }
                        if (!is_neighbor_inside_bb){
                            ComputeSeparatingFaces(curr_vox_tuple,q, face_1, face_2, false);
                            outer_faces.push_back(face_1);
                            outer_faces.push_back(face_2);
                        }
                    }
                }
            }
        }
    }
    
    num_faces_cage_ = outer_faces.size();
    cout << "Number of faces after cage extraction : " << num_faces_cage_ << endl;
    F_cage_ = MatrixXi::Zero(num_faces_cage_,3);
    
    // Construct mesh based on outer faces
    vector<Vector3d> vector_vertices;
    vector_vertices.clear();
    VectorXi added_vertices = - VectorXi::Ones(grid_face_vertices_.rows());
    int count = 0;
    for (int f = 0; f < num_faces_cage_; f++){
        Vector3i face = outer_faces[f];
        Vector3i new_face = Vector3i::Zero(3);
        for (int i = 0; i < 3; i++){
            int ind = face(i);
            int curr_vert_ind = added_vertices(ind);
            if (curr_vert_ind == -1){
                added_vertices(ind) = count;
                new_face(i) = count;
                vector_vertices.push_back(grid_face_vertices_.row(ind));
                count++;
            }
            else{
                new_face(i) = curr_vert_ind;
            }
        }
        F_cage_.row(f) = new_face;
    }
    num_vertices_cage_ = count;
    V_cage_ = MatrixXd::Zero(num_vertices_cage_,3);
    for (int i = 0; i < num_vertices_cage_; i++){
        V_cage_.row(i) = vector_vertices.at(i);
    }
}

void CageGenerator::FillingAlgorithm(const Vector3i &curr_vox_tuple, vector<Vector3i> &outer_faces) {
    Vector3i neighbor_vox_tuple;
    int i = curr_vox_tuple(0), j = curr_vox_tuple(1), k = curr_vox_tuple(2);
    int i_neighbor, j_neighbor, k_neighbor;
    int neighbor_partition_info;
    bool is_neighbor_inside_bb;
    Vector3i face_1, face_2;
    // Mark as visited
    voxel_partition_info_[i][j][k] = 1;
    for (int q = 0; q < NUM_FACES_CUBE; q++){
        neighbor_vox_tuple = curr_vox_tuple + Vector3i(voxel_face_normals_.row(q).cast<int>());
        // Check that neighboring voxel is inside the BB
        is_neighbor_inside_bb = true;
        for (int l = 0; l < 3; l++){
            if (neighbor_vox_tuple(l) < 0 || neighbor_vox_tuple(l) > n_vox[l] - 1){
                is_neighbor_inside_bb = false;
                break;
            }
        }
        if (is_neighbor_inside_bb){
            i_neighbor = neighbor_vox_tuple(0);
            j_neighbor = neighbor_vox_tuple(1);
            k_neighbor = neighbor_vox_tuple(2);
            neighbor_partition_info = voxel_partition_info_[i_neighbor][j_neighbor][k_neighbor];
            if (neighbor_partition_info == 1){
                // already visited
                continue;
            }
            else if(neighbor_partition_info == 0){
                // Feature voxel : Add the separating face
                ComputeSeparatingFaces(curr_vox_tuple,q, face_1, face_2, true);
                outer_faces.push_back(face_1);
                outer_faces.push_back(face_2);
            }
            else{
                // Visit neighbor
                FillingAlgorithm(neighbor_vox_tuple, outer_faces);
            }
        }
    }
}

void CageGenerator::SimplifyAdjacentFaces() {
    // Compute vertex-face adjacency vector
    VectorXi deleted_vertices = -VectorXi::Ones(num_vertices_cage_);
    VectorXi deleted_faces = -VectorXi::Ones(num_faces_cage_);
    vector < vector <int> > vertex_face_adjacency;
    vector <int> adj_faces;
    vector < Vector3i > new_faces;
    vector < int > corners;
    Vector3d current_point, adjface_point;
    int plane_dim;
    int other_dims[2] = {-1,-1};
    bool same_plane, deleted;
    for (int i = 0; i < num_vertices_cage_; i++) {
        vertex_face_adjacency.push_back(vector <int>());
    }
    for (int f = 0; f < num_faces_cage_; f++) {
        for (int j = 0; j < 3; j++){
            vertex_face_adjacency[F_cage_(f,j) ].push_back(f);
        }
    }
    for (int i = 0; i < num_vertices_cage_; i++) {
        current_point = V_cage_.row(i);
        adj_faces.clear();
        adj_faces = vertex_face_adjacency[i];
        deleted = false;
        if (adj_faces.size() == 6) {
            plane_dim = -1;
            // Find a dimension where the 6 neighbors lie in the same plane
            for (int d = 0; d < 3; d++){
                same_plane = true;
                // test each face
                for (int j = 0; j < 6; j++) {
                    // test each vertex of each face
                    for (int k = 0; k < 3; k++) {
                        adjface_point = V_cage_.row(F_cage_(adj_faces[j],k));
                        if ( abs(adjface_point(d) - current_point(d)) > epsilon_tests[d] ) {
                            same_plane = false;
                            break;
                        }
                    }
                    if (same_plane == false) break;
                }
                
                if (same_plane) {
                    plane_dim = d;
                    break;
                }
            }
            if (plane_dim != -1) {
                int count = 0;
                for (int k = 0 ; k < 3; k++){
                    if (k != plane_dim) {
                        other_dims[count] = k;
                        count++;
                    }
                }
                // DELETE the current vertex : retriangulate
                vector<int> neighbor_ids;
                deleted = true;
                int neighbor_id;
                for (int q = 0; q < 6; q++) {
                    int f = adj_faces[q];
                    for (int j = 0; j < 3; j++) {
                        neighbor_id = F_cage_(f,j);
                        if (neighbor_id != i && find(neighbor_ids.begin(), neighbor_ids.end(), neighbor_id) == neighbor_ids.end()) {
                            if (deleted_vertices(neighbor_id) == 1) {
                                deleted = false;
                            }
                            neighbor_ids.push_back(neighbor_id);
                        }
                    }
                    if (!deleted) break;
                }
                
                if (deleted) {
                    // Retriangulate : add new faces based on type of neighboring vertex
                    vector<int> ordered_neighbors;
                    int d0 = other_dims[0], d1 = other_dims[1];
                    for (int q = 0; q < 6; q++) {
                        if ( abs(V_cage_(neighbor_ids[q],d0) - current_point(d0)) > epsilon_tests[d0] &&
                            abs(V_cage_(neighbor_ids[q],d1) - current_point(d1)) < epsilon_tests[d1]) {
                            ordered_neighbors.push_back(neighbor_ids[q]);
                        }
                    }
                    for (int q = 0; q < 6; q++) {
                        if ( abs(V_cage_(neighbor_ids[q],d0) - current_point(d0)) < epsilon_tests[d0] &&
                            abs(V_cage_(neighbor_ids[q],d1) - current_point(d1)) > epsilon_tests[d1]) {
                            ordered_neighbors.push_back(neighbor_ids[q]);
                        }
                    }
                    for (int q = 0; q < 6; q++) {
                        if ( abs(V_cage_(neighbor_ids[q],d0) - current_point(d0)) > epsilon_tests[d0] &&
                            abs(V_cage_(neighbor_ids[q],d1) - current_point(d1)) > epsilon_tests[d1]) {
                            ordered_neighbors.push_back(neighbor_ids[q]);
                        }
                    }
                    if (ordered_neighbors.size() == 6) {
                        Vector3i face1(ordered_neighbors[0],ordered_neighbors[1],ordered_neighbors[2]);
                        Vector3i face2(ordered_neighbors[0],ordered_neighbors[1],ordered_neighbors[3]);
                        Vector2i adj5 = Vector2i::Zero(2), adj4 = Vector2i::Zero(2);
                        Vector3d v5 = V_cage_.row(ordered_neighbors[5]);
                        Vector3d vAdj;
                        int c = 0;
                        for (int q = 0; q < 5; q++){
                            vAdj = V_cage_.row(ordered_neighbors[q]);
                            if (abs(vAdj(d0) - v5(d0)) < res_vox[d0] + epsilon_tests[d0]  &&
                                abs(vAdj(d1) - v5(d1)) < res_vox[d1] + epsilon_tests[d1] ) {
                                adj5(c) = ordered_neighbors[q];
                                c++;
                            }
                        }
                        if (c != 2) continue;
                        c = 0;
                        for (int q = 0; q < 4; q++){
                            if (ordered_neighbors[q] != adj5(0) && ordered_neighbors[q] != adj5(1)) {
                                adj4(c) = ordered_neighbors[q];
                                c++;
                            }
                        }
                        Vector3i face3(ordered_neighbors[5],adj5(0), adj5(1));
                        Vector3i face4(ordered_neighbors[4],adj4(0), adj4(1));
                        new_faces.push_back(face1);
                        new_faces.push_back(face2);
                        new_faces.push_back(face3);
                        new_faces.push_back(face4);
                        
                        for (int f = 0; f < 6; f++) {
                            deleted_faces(adj_faces[f]) = 1;
                        }
                        deleted_vertices(i) = 1;
                    }
                }
            }
        }

    }
    // add faces that haven't been deleted
    for (int f = 0; f < num_faces_cage_; f++) {
        if (deleted_faces(f) != 1) {
            new_faces.push_back(F_cage_.row(f));
        }
    }
    
    num_faces_cage_ = new_faces.size();
    F_cage_ = MatrixXi::Zero(num_faces_cage_,3);
    cout << "Number of faces after mesh simplification : " << num_faces_cage_ << endl;
    // Construct mesh based on outer faces
    vector<Vector3d> vector_vertices;
    vector_vertices.clear();
    VectorXi added_vertices = - VectorXi::Ones(num_vertices_cage_);
    int count = 0;
    for (int f = 0; f < num_faces_cage_; f++){
        Vector3i face = new_faces[f];
        Vector3i new_face = Vector3i::Zero(3);
        for (int i = 0; i < 3; i++){
            int ind = face(i);
            int curr_vert_ind = added_vertices(ind);
            if (curr_vert_ind == -1){
                added_vertices(ind) = count;
                new_face(i) = count;
                vector_vertices.push_back(V_cage_.row(ind));
                count++;
            }
            else{
                new_face(i) = curr_vert_ind;
            }
        }
        F_cage_.row(f) = new_face;
    }
    num_vertices_cage_ = count;
    V_cage_ = MatrixXd::Zero(num_vertices_cage_,3);
    for (int i = 0; i < num_vertices_cage_; i++){
        V_cage_.row(i) = vector_vertices.at(i);
    }
}

void CageGenerator::ComputeCage(MatrixXd& bb_vertices, MatrixXi& bb_faces,
                                vector<MatrixXd>& feat_voxels_vertices, vector<MatrixXi>& feat_voxels_faces) {
    ComputePrincipalDirectionsBoundingBox(bb_vertices, bb_faces);
    InitializeVoxelPositions();
    FindFeatureVoxels(feat_voxels_vertices, feat_voxels_faces);
    ExtractOuterSurface();
    //SimplifyAdjacentFaces();
    SmoothCage();
}

int CageGenerator::TupleToIndex(const int& i, const int& j, const int& k) {
    int res = k * (n_vox[0]+1) * (n_vox[1]+1) + j * (n_vox[0]+1) + i;
    return res;
}

MatrixXd CageGenerator::GetGridVertices() {
    MatrixXd res = (pca_basis_matrix_ * grid_face_vertices_.transpose()).transpose();
    for (int i = 0; i < res.rows(); i++){
        res.row(i) = barycenter_ + Vector3d(res.row(i));
    }
    return res;
}

MatrixXi CageGenerator::GetCageFaces() { 
    return F_cage_;
}

void CageGenerator::SetSmoothingParameters(const int& n_i, const float& lambda) {
    lambda_smooth_ = lambda;
    num_smoothing_iterations_ = n_i;
}

void CageGenerator::SmoothCage() { 
    igl::cotmatrix(V_cage_,F_cage_,laplace_beltrami_matrix_);
    MatrixXd V_cage_prev = V_cage_;
    V_cage_smooth_ = V_cage_;
    SparseMatrix<double> S(num_vertices_cage_,num_vertices_cage_);
    S.setIdentity();
    BiCGSTAB< SparseMatrix<double> > solverImpl;
    if (implicit_smoothing){
        S = S - lambda_smooth_ * laplace_beltrami_matrix_;
        solverImpl.compute(S);
    }
    else S = S + lambda_smooth_ * laplace_beltrami_matrix_;
    Vector3d s, dir, dir_norm;
    vector<igl::Hit> hits;
    const double delta_T = 0.5;
    for (int i = 0; i < num_smoothing_iterations_; i++){
        if (implicit_smoothing)  V_cage_smooth_ =  solverImpl.solve(V_cage_prev);
        else V_cage_smooth_ = S * V_cage_prev;
        // Test for intersection with the mesh
        
        for (int j = 0; j < num_vertices_cage_; j++) {
            hits.clear();
            s = V_cage_prev.row(j);
            dir = V_cage_smooth_.row(j) - V_cage_prev.row(j);
            dir_norm = dir.normalized();
            if (ray_mesh_intersect(s,dir_norm,V_pca_basis_,F_mesh_,hits)){
                if (hits[0].t < dir.norm()) {
                V_cage_smooth_.row(j) = Vector3d(V_cage_prev.row(j)) + delta_T * hits[0].t * dir_norm;
                }
            }
        }
        V_cage_prev = V_cage_smooth_;
    }
}

MatrixXd CageGenerator::GetCage() {
    MatrixXd result = (pca_basis_matrix_ * V_cage_.transpose()).transpose();
    for (int i = 0 ; i < result.rows(); i++){
        result.row(i) = barycenter_ + Vector3d(result.row(i));
    }
    return result;
}

MatrixXd CageGenerator::GetSmoothCage() {
    MatrixXd result = (pca_basis_matrix_ * V_cage_smooth_.transpose()).transpose();
    for (int i = 0 ; i < result.rows(); i++){
        result.row(i) = barycenter_ + Vector3d(result.row(i));
    }
    return result;
}

















