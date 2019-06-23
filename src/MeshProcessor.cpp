//
//  MeshProcessor.cpp
//  3DA_project_CageBasedDef_bin
//
//  Created by Benjamin Barral on 12/05/2019.
//

#include "MeshProcessor.hpp"

// A class that computes the barycenter and extreme coordinates of a mesh

MeshProcessor::MeshProcessor(const MatrixXd &vertices, const MatrixXi &faces) { 
    V_ = vertices;
    F_ = faces;
    num_vertices_ = vertices.rows();
    num_faces_ = faces.rows();
    
    barycenter_ = Vector3d::Zero(3);
    for (int i = 0; i < num_vertices_; i++){
        barycenter_ += vertices.row(i) / num_vertices_;
    }
    max_vertices_ = new ExtremeVertex[3];
    min_vertices_ = new ExtremeVertex[3];
    ComputeExtremeVertices();
}

Vector3d MeshProcessor::GetBarycenter() {
    return barycenter_;
}

void MeshProcessor::ComputeExtremeVertices() {
    // Compute the index of the extreme vertices in all directions
    double minCoords[3] = {__DBL_MAX__, __DBL_MAX__, __DBL_MAX__};
    double maxCoords[3] = {__DBL_MIN__, __DBL_MIN__, __DBL_MIN__};
    Vector3d point;
    for (int i = 0; i < num_vertices_; i++){
        point = V_.row(i);
        for (int j = 0; j < 3; j++){
            if (point(j) > maxCoords[j]){
                maxCoords[j] = point(j);
                max_vertices_[j].vertex_index_ = i;
            }
            if (point(j) < minCoords[j]){
                minCoords[j] = point(j);
                min_vertices_[j].vertex_index_ = i;
            }
        }
    }
    for (int j = 0; j < 3; j++){
        max_vertices_[j].coord_ = maxCoords[j];
        min_vertices_[j].coord_ = minCoords[j];
    }
}

int MeshProcessor::GetMaximumVertexIndex(const int &dim) {
    return max_vertices_[dim].vertex_index_;
}

int MeshProcessor::GetMinimumVertexIndex(const int &dim) {
    return min_vertices_[dim].vertex_index_;
}

double MeshProcessor::GetBoundingBoxSize() { 
    double res = 0.;
    for (int j = 0; j < 3; j++){
        res += pow(max_vertices_[j].coord_ - min_vertices_[j].coord_,2);
    }
    return sqrt(res);
}





