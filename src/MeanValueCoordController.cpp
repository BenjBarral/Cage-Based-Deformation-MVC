//
//  MeanValueCoordController.cpp
//  3DA_project_CageBasedDef_bin
//
//  Created by Benjamin Barral on 11/05/2019.
//

#include "MeanValueCoordController.hpp"

static const double EPSILON_PLANAR_RATIO = 0.000001;

MeanValueCoordController::MeanValueCoordController() {
}

MeanValueCoordController::MeanValueCoordController(const MatrixXd &vertices_mesh, const MatrixXd &vertices_cage,
                                                   const MatrixXi &faces_mesh, const MatrixXi &faces_cage,
                                                   const double& bb_size) {
    V_mesh = vertices_mesh;
    V_cage = vertices_cage;
    V_cage_deformed = vertices_cage;
    F_mesh = faces_mesh;
    F_cage = faces_cage;
    num_vertices_mesh = V_mesh.rows();
    num_vertices_cage = V_cage.rows();
    num_faces_mesh = F_mesh.rows();
    num_faces_cage = F_cage.rows();
    epsilon_ = EPSILON_PLANAR_RATIO * bb_size;
}

void MeanValueCoordController::ComputeMVWeights() {
    mV_weights = MatrixXd::Zero(num_vertices_mesh, num_vertices_cage);
    VectorXd vector_D;
    MatrixXd matrix_U;
    Vector3d x, pJ, p1, p2, p3;
    double dJ;
    Vector3i indices;
    int num_planaer_cases1 = 0, num_planaer_cases2 = 0;
    for (int k = 0; k < num_vertices_mesh; k++){
        VectorXd weightsK = VectorXd::Zero(num_vertices_cage);
        double totalW = 0.;
        x = V_mesh.row(k);
        matrix_U = MatrixXd::Zero(num_vertices_cage,3);
        vector_D = VectorXd::Zero(num_vertices_cage);
        for (int j = 0; j  < num_vertices_cage; j++){
            pJ = V_cage.row(j);
            dJ = (pJ - x).norm();
            vector_D(j) = dJ;
            matrix_U.row(j) = (pJ - x) / dJ;
        }
        for (int f = 0; f < num_faces_cage; f++){
            indices = Vector3i(F_cage.row(f));
            p1 = V_cage.row(indices(0));
            p2 = V_cage.row(indices(1));
            p3 = V_cage.row(indices(2));
            double thetaIs[3] = {0.,0.,0.};
            double cIs[3] = {0.,0.,0.};
            double sIs[3] = {0.,0.,0.};
            double h = 0.;
            MatrixXd u123 = MatrixXd::Zero(3,3);
            for (int i = 0; i < 3; i++){
                int iMinus = (i-1) % 3;
                iMinus = (iMinus < 0) ? iMinus+3 : iMinus;
                Vector3d uIPlus = matrix_U.row(indices((i+1) % 3));
                Vector3d uIMinus = matrix_U.row(indices(iMinus));
                double lI = (uIPlus - uIMinus).norm();
                thetaIs[i] = 2. * asin(lI / 2.);
                h += thetaIs[i]/2.;
                u123.row(i) = matrix_U.row(indices(i));
            }
            if (M_PI - h < epsilon_){
                weightsK = VectorXd::Zero(num_vertices_cage);
                totalW = 0.;
                for (int i = 0; i < 3; i++){
                    int iMinus = (i-1) % 3;
                    iMinus = (iMinus < 0) ? iMinus+3 : iMinus;
                    //sin[θi]di−1di+1
                    double wI = sin(thetaIs[i]) * vector_D(indices((i+1)%3)) * vector_D(indices(iMinus));
                    if (isnan(wI)){
                        cout << "NaN" << endl;
                    }
                    weightsK(indices(i)) = wI;
                    totalW += wI;
                    
                }
                num_planaer_cases1 += 1;
                break;
            }
            double signDet = u123.determinant();
            if (signDet == 0){
                cout << "u123 = " << endl;
                cout << "(" << u123(0,0) << ", " << u123(0,1) << ", " << u123(0,2) << ")" << endl;
                cout << "(" << u123(1,0) << ", " << u123(1,1) << ", " << u123(1,2) << ")" << endl;
                cout << "(" << u123(2,0) << ", " << u123(2,1) << ", " << u123(1,2) << ")" << endl;
            }
            signDet = signDet / abs(signDet);
            bool discardTriangle = false;
            for (int i = 0; i < 3; i++){
                int iMinus = (i-1) % 3;
                iMinus = (iMinus < 0) ? iMinus+3 : iMinus;
                double cI = -1. + 2. * sin(h) * sin(h-thetaIs[i]) /
                (sin (thetaIs[(i+1)%3]) * sin(thetaIs[iMinus]));
                cIs[i] = cI;
                if (cI < -1.) {
                    cI = -1.;
                    cout << "cI < -1 " << endl;
                }
                if (cI > 1.){
                    cI = 1.;
                    cout << "cI > 1 " << endl;
                }
                double sI = signDet * sqrt(1. - cI*cI);
                if (isnan(sI)) {
                    cout << "NaN" << endl;
                }
                if (abs(sI) < epsilon_){
                    discardTriangle = true;
                    num_planaer_cases2 +=1;
                    break;
                }
                else sIs[i] = sI;
            }
            if (!discardTriangle){
                for (int i = 0; i < 3; i++){
                    int iPlus = (i+1) % 3;
                    int iMinus = (i-1) % 3;
                    iMinus = (iMinus < 0) ? iMinus+3 : iMinus;
                    double dI = vector_D(indices(i));
                    double wI = (thetaIs[i] - cIs[iPlus] * thetaIs[iMinus] - cIs[iMinus] * thetaIs[iPlus] ) /
                    (dI * sin(thetaIs[iPlus]) * sIs[iMinus]);
                    if (isnan(wI)){
                        cout << "NaN" << endl;
                    }
                    weightsK(indices(i)) += wI;
                    totalW += wI;
                }
            }
        }
        weightsK /= totalW;
        mV_weights.row(k) = weightsK;
    }
    cout << "[Mean Value Coordinates] number of planar cases : " << num_planaer_cases1 << " and " << num_planaer_cases2 << endl;
}

void MeanValueCoordController::SetDeformedCage(const MatrixXd& vertices_cage_deformed){
    V_cage_deformed = vertices_cage_deformed;
}

MatrixXd MeanValueCoordController::MVInterpolate(){
    return mV_weights * V_cage_deformed;
}



