//
//  main.cpp
//  3DA_project_CageBasedDef_bin
//
//  Created by Benjamin Barral on 05/02/2019.
//

#include <stdio.h>
#define DIM_X 0
#define DIM_Y 1
#define DIM_Z 2


#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <iostream>
#include <igl/readPLY.h>
#include <igl/writePLY.h>
#include <igl/writeOBJ.h>
#include <igl/file_exists.h>
#include <Eigen/Geometry>
#include <cmath>
#include <iterator>
#include "MeanValueCoordController.hpp"
#include "MeshProcessor.hpp"
#include <random>
#include "CageGenerator.hpp"
#include <igl/opengl/glfw/Viewer.h>
#include "DeformCageViewerPlugin.hpp"


using namespace std;
using namespace Eigen;
using namespace igl::opengl::glfw;


string mesh_file_name = "../../ HumanBody_3670.ply";
string cage_file_name = "../../Human_Cage_V5.obj";

const bool compute_automatic_cage = false;
float sparseness_cage = 0.5; // For automatic cage generation : CHANGE THIS for a sparser or denser cage : the larger the parameter the denser the cage

int main(int argc, char *argv[])
{
    srand(time(NULL));
    clock_t start;
    
    // LOAD MESHES
    MatrixXd V_mesh,V_cage;
    MatrixXi F_mesh,F_cage;
    if (mesh_file_name.find("obj") != string::npos) igl::readOBJ(mesh_file_name,V_mesh,F_mesh);
    else if (mesh_file_name.find("ply") != string::npos) igl::readPLY(mesh_file_name,V_mesh,F_mesh);
    else if (mesh_file_name.find("off") != string::npos) igl::readOFF(mesh_file_name,V_mesh,F_mesh);
    else {
        cout << "Mesh file not recognized " << endl;
        return;
    }
    
    // Init the viewer
    Viewer viewer;
    viewer.core.is_animating = true;
    viewer.append_mesh();
    viewer.append_mesh();
    float point_size = 18.;
    
    // Generate automatic cage
    MatrixXd V_cage_automatic, V_cage_automatic_smooth;
    float lambda_smooth_implicit = .6;
    int num_iterations_smoothing = 2;
    CageGenerator cage_generator(V_mesh, F_mesh, num_iterations_smoothing, lambda_smooth_implicit, sparseness_cage);
    if (compute_automatic_cage) {
        // Generate cage
        MatrixXd bb_vertices;
        MatrixXi bb_faces;
        vector<MatrixXd> feat_voxels_vertices;
        vector<MatrixXi> feat_voxels_faces;
        start = clock();
        cage_generator.ComputeCage(bb_vertices, bb_faces, feat_voxels_vertices, feat_voxels_faces);
        cout << "Time to compute cage : " << (clock() - start) / (double) CLOCKS_PER_SEC << endl;
        
        V_cage_automatic = cage_generator.GetCage();
        V_cage_automatic_smooth = cage_generator.GetSmoothCage();
        //MatrixXd V_cage_automatic = cage_generator.GetSmoothCage();
        MatrixXi F_cage_automatic = cage_generator.GetCageFaces();
        // Show bounding box
        viewer.append_mesh();
        viewer.data_list[3].set_mesh(bb_vertices, bb_faces);
         viewer.data_list[3].show_faces = false;
        
        // Show automatically generated cage
        V_cage = V_cage_automatic_smooth;
        F_cage = F_cage_automatic;
        
        // Save mesh in PLY file
        igl::writePLY("../../GeneratedCage.ply", V_cage_automatic_smooth, F_cage_automatic);
    }
    else {
        if (cage_file_name.find("obj") != string::npos) igl::readOBJ(cage_file_name,V_cage,F_cage);
        else if (cage_file_name.find("ply") != string::npos) igl::readPLY(cage_file_name,V_cage,F_cage);
        else if (cage_file_name.find("off") != string::npos) igl::readOFF(cage_file_name,V_cage,F_cage);
        else {
            cout << "Cage file not recognized " << endl;
            return;
        }
    }
    
    // Get barycenter and extreme points of mesh
    MeshProcessor cage_processor(V_cage,F_cage);
    Vector3d cage_barycenter = cage_processor.GetBarycenter();
    MeshProcessor mesh_processor(V_mesh,F_mesh);
    Vector3d mesh_barycenter = mesh_processor.GetBarycenter();
    int dim = DIM_X;
    int max_vertY_ind = cage_processor.GetMaximumVertexIndex(dim);
    double bb_size = mesh_processor.GetBoundingBoxSize();
    
    // Print out number of triangles
    int num_faces_mesh = F_mesh.rows();
    int num_vertices_cage = V_cage.rows();
    int num_faces_cage = F_cage.rows();
    int num_vertices_mesh = V_mesh.rows();
    cout << "Cage : " << num_faces_cage << " triangles, " << num_vertices_cage << " vertices." << endl;
    cout << "Mesh : " << num_faces_mesh << " triangles, " << num_vertices_mesh << " vertices." << endl;
    
    
    // Visualize axes
    MatrixXd axes_points1 = MatrixXd::Zero(3, 3);
    MatrixXd axes_points2 = MatrixXd::Zero(3, 3);
    axes_points2.coeffRef(0, 0) = 1;
    axes_points2.coeffRef(1,1) = 1;
    axes_points2.coeffRef(2,2) = 1;
    MatrixXd axes_colors = axes_points2;
    axes_points2 *= 0.25;
    double offset_z = 1.3;
    for (int i = 0; i<3; i++){
        axes_points2.coeffRef(i, 2) += offset_z;
        axes_points1.coeffRef(i, 2) += offset_z;
    }
    viewer.data_list[0].add_edges(axes_points1, axes_points2, axes_colors);
    
    viewer.data_list[1].set_mesh(V_mesh, F_mesh);
    viewer.data_list[1].show_lines = false;
    MatrixXd cage_points_colors = MatrixXd::Ones(num_vertices_cage,3);
    viewer.data_list[2].set_mesh(V_cage, F_cage);
    viewer.data_list[2].show_faces = false;
    viewer.data_list[2].add_points(V_cage, cage_points_colors);
    viewer.data_list[2].point_size = point_size;
            
     
     // Mean Value Coordinates
     MeanValueCoordController mVCoord_controller(V_mesh,V_cage,F_mesh,F_cage, bb_size);
     start = clock();
     //if (!compute_automatic_cage)
         mVCoord_controller.ComputeMVWeights();
     cout << "Time to compute weights : " << (clock() - start) / (double) CLOCKS_PER_SEC << endl;
     start = clock();
    MatrixXd V_mesh_deformed = V_mesh;
     //if (!compute_automatic_cage)
         V_mesh_deformed = mVCoord_controller.MVInterpolate();
     float timeInterp = (clock() - start) / (double) CLOCKS_PER_SEC;
     float fps = 1. / timeInterp;
     cout << "Time to interpolate : " << timeInterp << endl;
     cout << "fps = " << fps << endl;
     
     MatrixXd V_cage_deformed = V_cage;
     mVCoord_controller.SetDeformedCage(V_cage_deformed);
    
    
    // Attach a plugin to handle deformation interaction
    int cage_data_id = 2, mesh_data_id = 1;
    DeformCageViewerPlugin def_cage_plugin(V_cage, cage_points_colors, F_cage,
                                           mesh_data_id, cage_data_id, mVCoord_controller);
    viewer.plugins.push_back(&def_cage_plugin);
    
    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    def_cage_plugin.init(&viewer);
    viewer.plugins.push_back(&menu);
    
    // Wave deformation
    float wave_time_start;
    bool wave_isActive = false;
    float wave_duration = .75;
    int wave_index_point;
    
    
    // Deform cage
    float cage_resize_ratioZ = 1.,  cage_resize_ratioY = 1., cage_resize_ratioX = 1.;
    int smooth_cage_slider = 1;
    
    // UI instructions
    cout << endl;
    cout << "CAGE INTERACTION : " << endl;
    cout << "J + click on a face : select/unselect a joint point (appears in blue)" << endl;
    cout << "S + click on a face : " << endl;
    cout << "   If joint point selected : select group of points at distance lesser than joint point from the click (appear in green) by rotating around the joint point" << endl;
    cout << "   If no joint point selected : move one point (appears in red) " << endl;
    cout << "P : Hide points" << endl;
    
    menu.callback_draw_custom_window = [&]()
    {
        bool update_cage = false;
        bool update_mesh = false;
        // Define next window position + size
        ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(800, 100), ImGuiSetCond_FirstUseEver);
        ImGui::Begin(
                     "MyProperties", nullptr,
                     ImGuiWindowFlags_NoSavedSettings
                     );
        ImGui::End();
        
        
        if (ImGui::Button("Show interpolated")){
            viewer.data_list[1].clear();
            viewer.data_list[1].set_mesh(V_mesh_deformed, F_mesh);
            viewer.data_list[2].clear();
            viewer.data_list[2].set_mesh(V_cage_deformed, F_cage);
        }
        if (ImGui::Button("Show original")){
            viewer.data_list[1].clear();
            viewer.data_list[1].set_mesh(V_mesh, F_mesh);
            viewer.data_list[2].clear();
            viewer.data_list[2].set_mesh(V_cage, F_cage);
        }
        if (ImGui::SliderFloat("Stretching cage on Z axis", &cage_resize_ratioZ, 0, 2)){
            for (int i = 0; i < num_vertices_cage; i++){
                V_cage_deformed(i,DIM_Z) = cage_barycenter(DIM_Z) +
                (V_cage(i,DIM_Z) -  cage_barycenter(DIM_Z)) * cage_resize_ratioZ;
            }
            update_mesh = true;
            update_cage = true;
        }
        if (ImGui::SliderFloat("Stretching cage on Y axis", &cage_resize_ratioY, 0, 2)){
            for (int i = 0; i < num_vertices_cage; i++){
                V_cage_deformed(i,DIM_Y) = cage_barycenter(DIM_Y) +
                (V_cage(i,DIM_Y) -  cage_barycenter(DIM_Y)) * cage_resize_ratioY;
            }
            update_mesh = true;
            update_cage = true;
        }
        if (ImGui::SliderFloat("Stretching cage on X axis", &cage_resize_ratioX, 0, 2)){
            for (int i = 0; i < num_vertices_cage; i++){
                V_cage_deformed(i,DIM_X) = cage_barycenter(DIM_X) +
                (V_cage(i,DIM_X) -  cage_barycenter(DIM_X)) * cage_resize_ratioX;
            }
            update_mesh = true;
            update_cage = true;
        }
     
        if (compute_automatic_cage) {
            if (ImGui::SliderInt("Show coarse/smooth cage", &smooth_cage_slider, 0, 1)){
                V_cage_deformed = (smooth_cage_slider == 1) ? V_cage_automatic_smooth : V_cage_automatic;
                update_cage = true;
            }
            ImGui::SliderFloat("lamba smoothing", &lambda_smooth_implicit, 0,2);
            ImGui::SliderInt("Smoothing iterations", &num_iterations_smoothing, 1, 15);
            if (ImGui::Button("Smooth Cage")) {
                cage_generator.SetSmoothingParameters(num_iterations_smoothing, lambda_smooth_implicit);
                cage_generator.SmoothCage();
                V_cage_automatic_smooth = cage_generator.GetSmoothCage();
                V_cage_deformed = V_cage_automatic_smooth;
                update_cage = true;
            }
        }
        
        if (ImGui::Button("Reset")){
            V_cage_deformed = V_cage;
            V_mesh_deformed = V_mesh;
            viewer.data_list[1].set_vertices(V_mesh_deformed);
            update_cage = true;
        }
        
        if(ImGui::Button("Random point wave deformation")){
            wave_index_point = rand() % num_vertices_cage;
            wave_isActive = true;
            wave_time_start = clock();
        }
        ImGui::SliderFloat("Wave duration", &wave_duration, 0.2, 2);
        
        if (wave_isActive){
            
            float t = 2. * (clock() - wave_time_start) /((double) CLOCKS_PER_SEC * wave_duration);
            if (t > 1.){
                wave_isActive = false;
                V_cage_deformed.row(wave_index_point) = V_cage.row(wave_index_point);
            }
            else{
                double ratio = 1. + sin(M_PI * t);
                Vector3d original_vert = V_cage.row(wave_index_point);
                Vector3d deformed_cage_vert =  cage_barycenter +
                ratio * (original_vert - cage_barycenter);
                V_cage_deformed.row(wave_index_point) = deformed_cage_vert;
            }
            update_mesh = true;
            update_cage = true;
        }
        
        if (update_cage){
            viewer.data_list[2].set_vertices(V_cage_deformed);
            viewer.data_list[2].set_points(V_cage_deformed, cage_points_colors);
            def_cage_plugin.ResetVCage(V_cage_deformed, cage_points_colors);
            mVCoord_controller.SetDeformedCage(V_cage_deformed);
        }
        if (update_mesh){
            V_mesh_deformed = mVCoord_controller.MVInterpolate();
            viewer.data_list[1].set_vertices(V_mesh_deformed);
        }
        
        
            
    };
    
    // Call GUI
    viewer.launch();
 
}


