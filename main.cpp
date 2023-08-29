#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include "happly.h"
#include "icp.h"
#include "utils.h"



int main(int argc, char* argv[]) {

    using namespace std;
    // load the ply file
    if(argc == 1){
        cout << "Kindly enter the location of the plyfile as the argument";
        return 1;
    }
    char* plyfile_location = argv[1];

    // get the plyfile location

    string plyfile = plyfile_location;// "/home/brazenparadox/Documents/MASc/thesis/code_repo/3d_registration_castellani_and_bartoli_python/demo/bunny/reconstruction/bun_zipper_res3.ply";

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

    Eigen::MatrixXf dataview(3, x.size());

    // rotate the pcl matrix by 40 degrees along z axis
    dataview = rotate_matrix_along_z(pclmatrix, 40);

    // rotate the matrix along x
    dataview = rotate_matrix_along_x(dataview, 23);

    //rotate the matrix along y
    //dataview = rotate_matrix_along_y(dataview, 0); -- no significance since the degree is zero. Can be uncommented if degree is non-zero.

    // translate the matrix by 2 and 3 in x and y axes respectively
    dataview = translate_matrix(dataview, 2, 3, 0);

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
        // cout << "\nThe singular values from svd are " << svd.singularValues();
        // cout << "\n The u matrix is" << svd.matrixU();
        // cout << "\n The v matrix is" << svd.matrixV();

        Eigen::Matrix<float, 3, 3> u = svd.matrixU();
        Eigen::Matrix<float, 3, 3> v = svd.matrixV();
        Eigen::DiagonalMatrix<float, 3, 3> s;

        int det_v = round(v.determinant());
        int det_u = round(u.determinant());

        int determinanant_prod = det_u * det_v;
        if (determinanant_prod == 1) {
            s.diagonal() << 1, 1, 1;
        } else if (determinanant_prod == -1) {
            //cout << "\n\nproduct of determinants is -1";
            s.diagonal() << 1, 1, -1;
        }

        // cout << "\n\n singular values are " << s.diagonal();

        // calcuation rotation matrix
        Eigen::Matrix<float, 3, 3> rotational_matrix;
        rotational_matrix = v * s * u.transpose();
        // cout << "\n\nrotational_matrix" << rotational_matrix;

        // calculate translation matrix
        Eigen::Matrix<float, 3, 1> translational_matrix = model_mean - rotational_matrix * data_mean;
        // cout << "\n\ntranslational_matrix" << translational_matrix;

        // transform matrix
        dataview = transform_matrix(dataview, rotational_matrix, translational_matrix);
        data_mean = calculate_centroid(dataview);
        sse = calculate_sse(dataview, pclmatrix);

        cout << "\nSSE is " << sse << " in iteration " << iteration++;
    }

}
