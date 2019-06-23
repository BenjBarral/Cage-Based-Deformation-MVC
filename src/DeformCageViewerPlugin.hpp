//
//  DeformCageViewerPlugin.hpp
//  3DA_project_CageBasedDef_bin
//
//  Created by Benjamin Barral on 17/05/2019.
//

// Class for handling user interaction with the cage

#ifndef DeformCageViewerPlugin_hpp
#define DeformCageViewerPlugin_hpp

#include <stdio.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include "MeanValueCoordController.hpp"

using namespace igl::opengl::glfw;
using namespace Eigen;

static const double MAX_EDGE_SUM_DISTANCES_RATIO_ = 3.;

class DeformCageViewerPlugin : public ViewerPlugin{
private:
    MatrixXd V_cage_, V_cage_deformed_, C_points_cage_;
    MatrixXi F_cage_;
    int num_vertices_, num_faces_;
    int viewer_data_cage_id_, viewer_data_mesh_id_;
    int current_selected_point_id_, current_joint_selected_id_;
    bool is_selecting_point_, is_moving_point_, is_selecting_joint_point_, joint_point_selected_;
    double z_selected_point_;
    MeanValueCoordController mVCoord_controller_;
    Vector3d joint_point_;
    vector < vector < int > > neighbors_;
    Vector3d initial_mouse_world_pos_;
    vector< int > points_to_move_;
    VectorXi visited_points_;
    
public:
    DeformCageViewerPlugin(){}
    DeformCageViewerPlugin(const MatrixXd& V, const MatrixXd& C,
                           const MatrixXi& F, const int& viewerMeshDataId, const int& viewerCageDataId,
                           const MeanValueCoordController& mVCoord_ctrlr) {
        V_cage_ = V;
        V_cage_deformed_ = V;
        F_cage_ = F;
        viewer_data_cage_id_ = viewerCageDataId;
        viewer_data_mesh_id_ = viewerMeshDataId;
        num_vertices_ = V.rows();
        num_faces_ = F.rows();
        C_points_cage_ = C;
        is_selecting_point_ = false;
        is_moving_point_ = false;
        joint_point_selected_ = false;
        is_selecting_joint_point_ = false;
        
        mVCoord_controller_ = mVCoord_ctrlr;
    
        for (int i = 0 ; i < num_vertices_; i++){
            neighbors_.push_back(vector < int >());
        }
        for (int f = 0; f < num_faces_; f++){
            int current_ind, neighbor_ind;
            for (int i = 0; i < 3; i++){
                current_ind = F_cage_(f,i);
                for (int j = 0; j < 3; j++){
                    if (j!=i){
                        neighbor_ind = F_cage_(f,j);
                        neighbors_[current_ind].push_back(neighbor_ind);
                    }
                }
            }
        }
        
        // TODO : print the key instructions
    }
    
    void ResetVCage(const MatrixXd& V, const MatrixXd& colors) {
        V_cage_ = V;
        V_cage_deformed_ = V;
        is_selecting_point_ = false;
        is_moving_point_ = false;
        is_selecting_joint_point_ = false;
        joint_point_selected_ = false;
        C_points_cage_ = colors;
    }
    
    bool mouse_down(int button, int modifier) {
        int face_id;
        Eigen::Vector3f barycentric_coord;
        // Ray casting
        double x = viewer->current_mouse_x;
        double y = viewer->core.viewport(3) - viewer->current_mouse_y;
        
        if (button == GLFW_MOUSE_BUTTON_LEFT && (is_selecting_point_ || is_selecting_joint_point_)) {
            if(igl::unproject_onto_mesh(Vector2f(x,y), viewer->core.view,
                                        viewer->core.proj, viewer->core.viewport, V_cage_deformed_, F_cage_, face_id, barycentric_coord))
            {
                double max = DBL_MIN;
                int i_max = -1;
                Vector3d hitPoint1 = Vector3d::Zero(3);
                for (int i = 0; i < 3; i++){
                    hitPoint1 += barycentric_coord(i) * V_cage_deformed_.row(F_cage_(face_id, i));
                    if (barycentric_coord(i) > max) {
                        max = barycentric_coord(i);
                        i_max = i;
                    }
                }
                int v_selected = F_cage_(face_id, i_max);
                if (is_selecting_joint_point_) {
                    if (joint_point_selected_) {
                        joint_point_selected_ = false;
                        C_points_cage_.row(current_joint_selected_id_) << 1,1,1;
                    }
                    else {
                        joint_point_ = V_cage_deformed_.row(v_selected); // CHANGE : ADD POINT 
                        C_points_cage_.row(v_selected) << 0,0,1;
                        current_joint_selected_id_ = v_selected;
                        joint_point_selected_ = true;
                    }
                    viewer->data_list[viewer_data_cage_id_].set_mesh(V_cage_deformed_, F_cage_);
                    viewer->data_list[viewer_data_cage_id_].set_points(V_cage_deformed_, C_points_cage_);
                }
                else {
                    if (joint_point_selected_) {
                        points_to_move_.clear();
                        visited_points_ = VectorXi::Zero(num_vertices_);
                        double max_dist = (Vector3d(V_cage_deformed_.row(v_selected)) - joint_point_).norm();
                        AddPointsToMove(v_selected,max_dist, v_selected, 0);
                        V_cage_ = V_cage_deformed_;
                        for (int i = 0; i < points_to_move_.size(); i++) {
                            int ind = points_to_move_[i];
                            C_points_cage_.row(ind) << 0,1,0;
                            viewer->data_list[viewer_data_cage_id_].set_points(V_cage_deformed_, C_points_cage_);
                        }
                    }
                    else {
                        current_selected_point_id_ = v_selected;
                        C_points_cage_.row(current_selected_point_id_) << 1,0,0;
                        viewer->data_list[viewer_data_cage_id_].set_points(V_cage_deformed_, C_points_cage_);
                    }
                    
                    double x = viewer->current_mouse_x;
                    double y = viewer->core.viewport(3) - viewer->current_mouse_y;
                    Vector2d mouse_pos = Vector2d(x,y);
                    Vector4d homo = joint_point_.homogeneous();
                    Matrix4d view_matrix = (viewer->core.view).cast<double>();
                    Vector4d view_coordinates_homo = view_matrix * homo;
                    Vector3d view_coordinates = view_coordinates_homo.hnormalized();
                    double z = view_coordinates(2);
                    z_selected_point_ = z;
                    initial_mouse_world_pos_ = MouseToWorld(mouse_pos);
                    is_moving_point_ = true;
                }
                return true;
            }
        }
        return false;
    }
    
    void AddPointsToMove(const int& current_ind, const double& max_dist,
                            const int& start_ind, const double & sum_dists){
        double dist = (V_cage_deformed_.row(current_ind) - V_cage_deformed_.row(start_ind)).norm();
        if (dist < max_dist && sum_dists < MAX_EDGE_SUM_DISTANCES_RATIO_ * max_dist){
            points_to_move_.push_back(current_ind);
            visited_points_(current_ind) = 1;
            vector <int> current_neighbors = neighbors_[current_ind];
            int neighbor_ind;
            for (int i = 0; i < current_neighbors.size(); i++){
                neighbor_ind = current_neighbors[i];
                if (visited_points_(neighbor_ind) == 0){
                    double new_sum_dists = sum_dists +
                    (V_cage_deformed_.row(current_ind) - V_cage_deformed_.row(neighbor_ind)).norm();
                    AddPointsToMove(neighbor_ind, max_dist, start_ind, new_sum_dists);
                }
            }
        }
    }
    
    Vector3d MouseToWorld(const Vector2d& mouse_pos) {
        Vector4d homo;
        // Ray casting
        Vector3d s,dir;
        igl::unproject_ray(mouse_pos, viewer->core.view, viewer->core.proj, viewer->core.viewport,s,dir);
        Vector4d dir_4d = (viewer->core.view).cast<double>() * dir.homogeneous();
        Vector3d dir_view = dir_4d.hnormalized(); // position in view world
        Vector3d current_view_position = dir_view * z_selected_point_ / dir_view(2);
        
        homo = ( (viewer->core.view).cast<double>() ).inverse() * current_view_position.homogeneous();
        return homo.hnormalized();
    }
    bool mouse_move(int button, int modifier) {
        if (is_moving_point_) {
            double x = viewer->current_mouse_x;
            double y = viewer->core.viewport(3) - viewer->current_mouse_y;
            Vector2d mouse_pos(x,y);
            Vector3d current_mouse_world_position = MouseToWorld(mouse_pos);
            Vector3d initial_v = (initial_mouse_world_pos_ - joint_point_).normalized();
            Vector3d new_v = (current_mouse_world_position - joint_point_).normalized();
            Vector3d new_w = (new_v.cross(initial_v)).normalized();
            Vector3d new_u = (new_v.cross(new_w)).normalized();
            Vector3d initial_u = (initial_v.cross(new_w)).normalized();
            
            if (joint_point_selected_){
                Vector3d new_position, initial_position;
                for (int i = 0; i < points_to_move_.size(); i++) {
                    int ind = points_to_move_[i];
                    initial_position = V_cage_.row(ind);
                    new_position = current_mouse_world_position +
                    (initial_position - initial_mouse_world_pos_).dot(initial_v) * new_v +
                    (initial_position - initial_mouse_world_pos_).dot(initial_u) * new_u +
                    (initial_position - initial_mouse_world_pos_).dot(new_w) * new_w;
                    V_cage_deformed_.row(ind) = new_position;
                }
            }
            else {
                V_cage_deformed_.row(current_selected_point_id_) = current_mouse_world_position;
            }
            viewer->data_list[viewer_data_cage_id_].set_mesh(V_cage_deformed_, F_cage_);
            viewer->data_list[viewer_data_cage_id_].set_points(V_cage_deformed_, C_points_cage_);
            MatrixXd P1(1,3), P2(1,3), C(1,3);
            P1.row(0) = joint_point_;
            P2.row(0) = current_mouse_world_position;
            C.row(0) << 1,0,0;
            //viewer->data_list[viewer_data_cage_id_].set_edges(P1, P2, C);
            
            mVCoord_controller_.SetDeformedCage(V_cage_deformed_);
            viewer->data_list[viewer_data_mesh_id_].set_vertices(mVCoord_controller_.MVInterpolate());
            return true;
        }
        return false;
    }
    
    bool mouse_up(int button, int modifier){
        if (button == GLFW_MOUSE_BUTTON_LEFT && is_selecting_point_) {
            if (joint_point_selected_) {
                for (int i = 0; i < points_to_move_.size(); i++) {
                    int ind = points_to_move_[i];
                    C_points_cage_.row(ind) << 1,1,1;
                }
            }
            else{
                C_points_cage_.row(current_selected_point_id_)<<1,1,1;
            }
            viewer->data_list[viewer_data_cage_id_].set_mesh(V_cage_deformed_, F_cage_);
            viewer->data_list[viewer_data_cage_id_].set_points(V_cage_deformed_, C_points_cage_);
            is_moving_point_ = false;
            return true;
        }
        return false;
    }
    
    bool key_down(int key, int modifiers) {
        if (key == 'S') {
            is_selecting_point_ = true;
            return true;
        }
        if (key == 'J') {
            is_selecting_joint_point_ = true;
            return true;
        }
        if (key == 'P') {
            viewer->data_list[viewer_data_cage_id_].clear();
            viewer->data_list[viewer_data_cage_id_].set_mesh(V_cage_deformed_, F_cage_);
            return true;
        }
        return false;
    }
    
    bool key_up(int key, int modifiers) {
        if (key == 'S') {
            is_selecting_point_ = false;
            return true;
        }
        if (key == 'J') {
            is_selecting_joint_point_ = false;
            return true;
        }
        return false;
    }
};

#endif /* DeformCageViewerPlugin_hpp */
