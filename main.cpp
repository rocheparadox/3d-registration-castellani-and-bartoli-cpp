#include <iostream>
#include <eigen3/Eigen/Dense>
#include "happly.h"
#include "icp.h"
#include <eigen3/Eigen/SVD>


int main() {

    using namespace std;
    // load the ply file
    string plyfile = "/home/brazenparadox/Documents/MASc/thesis/code_repo/3d_registration_castellani_and_bartoli_python/demo/bunny/reconstruction/bun_zipper_res3.ply";

    happly::PLYData plyIn(plyfile);
    std::vector<float> x = plyIn.getElement("vertex").getProperty<float>("x");
    std::vector<float> y = plyIn.getElement("vertex").getProperty<float>("y");
    std::vector<float> z = plyIn.getElement("vertex").getProperty<float>("z");

    // Create matrix from plyfile data
    Eigen::MatrixXf pclmatrix(3, x.size());
    for(int i=0; i<x.size(); i++){
        pclmatrix(0,i) = x[i];
    }
    for(int i=0; i<x.size(); i++){
        pclmatrix(1,i) = y[i];
    }
    for(int i=0; i<x.size(); i++){
        pclmatrix(2,i) = z[i];
    }
//    for(int i=0; i<x.size(); i++){
//        pclmatrix(3, i) = 1;
//    }

    Eigen::MatrixXf dataview(3, x.size());

    // rotate the matrix along z
    int degree = 40;
    float theta = degree * 3.14/180;
    Eigen::Matrix<float, 3, 3> zrotational_matrix;
    zrotational_matrix << cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1 ;
    dataview << zrotational_matrix*pclmatrix;

//    // rotate the matrix along x
//    degree = 23;
//    theta = degree * 3.14/180;
//    Eigen::Matrix<float, 3, 3> xrotational_matrix;
//    xrotational_matrix << 1, 0, 0, 0, cos(theta), -sin(theta), 0, sin(theta), cos(theta);
//    dataview << xrotational_matrix*dataview;
//
//    // rotate the matrix along y
//    degree = 13;
//    theta = degree * 3.14/180;
//    Eigen::Matrix<float, 3, 3> yrotational_matrix;
//    yrotational_matrix << cos(theta), 0, -sin(theta),  0, 1, 0, sin(theta), 0, cos(theta)  ;
//    dataview << yrotational_matrix*dataview;

    //cout << zrotational_matrix << endl;


    //cout << dataview << endl;

    // translate the matrix by 2 in x and 3 in y
//    Eigen::Matrix<float, 3, 3> translational_matrix;
//    translational_matrix << 1, 0, 0, 2,
//            0, 1, 0, 3,
//            0, 0, 1, 0,
//            0, 0, 0, 1;
//    dataview << translational_matrix*dataview;

    Eigen::Vector3f translation_matrix = Eigen::Vector3f(2.0, 3.0, 0.0);
    //cout << translation_matrix;
    dataview << dataview.colwise() + translation_matrix;

    // calculate centroid of model view
    Eigen::Matrix<float, 3, 1> model_mean = calculate_centroid(pclmatrix);
    //cout << model_mean;
    // calculate centroid of data view
    Eigen::Matrix<float, 3, 1> data_mean = calculate_centroid(dataview);

    float sse = 1.1; // arbitrary value initialization
    int iteration = 0;
    while (sse > 0.1) {
        // get correspondent points
        Eigen::MatrixXf correspondent_points(3, x.size());
        correspondent_points = get_correspondence_points(pclmatrix, dataview);
        // cout << "The correspondent points are " << correspondent_points;
        //calculate cross covariance matrix
        Eigen::Matrix<float, 3, 3> cross_covariance = calculate_cross_covariance(dataview, correspondent_points);

        // decompose using svd
        Eigen::BDCSVD<Eigen::Matrix<float, 3, 3>> svd(cross_covariance, Eigen::ComputeFullU | Eigen::ComputeFullV);
        // cout << "\ncross covariance is " << cross_covariance;
        //svd.computeU();
        //svd.computeV();
        // cout << "\nThe singular values are " << svd.singularValues();
        //cout << "\n The u matrix is" << svd.matrixU();
        //cout << "\n The v matrix is" << svd.matrixV();

        Eigen::Matrix<float, 3, 3> u = svd.matrixU();
        Eigen::Matrix<float, 3, 3> v = svd.matrixV();
        Eigen::DiagonalMatrix<float, 3, 3> s;

        int det_v = round(v.determinant());
        int det_u = round(u.determinant());

        int determinanant_prod = det_u * det_v;
        if (determinanant_prod == 1) {
            s.diagonal() << 1, 1, 1;
        } else if (determinanant_prod == -1) {
//            cout << "\n\nproduct of determinants is -1";
            s.diagonal() << 1, 1, -1;
        }

        // cout << "\n\n singular values are " << s.diagonal();

        // calcuation rotation matrix
        Eigen::Matrix<float, 3, 3> rotational_matrix;
        rotational_matrix = v * s * u.transpose();
        //cout << "\n\nrotational_matrix" << rotational_matrix;

        // calculate translation matrix

        Eigen::Matrix<float, 3, 1> translational_matrix = model_mean - rotational_matrix * data_mean;
        //cout << "\n\ntranslational_matrix" << translational_matrix;

        // transform matrix
        dataview = transform_matrix(dataview, rotational_matrix, translational_matrix);
        data_mean = calculate_centroid(dataview);
        sse = calculate_sse(dataview, pclmatrix);

        cout << "\nSSE is " << sse << " in iteration " << iteration++;
    }

}
