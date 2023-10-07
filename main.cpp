#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include "happly.h"
#include "icp.h"
#include "utils.h"

#include<chrono>

Eigen::MatrixXf get_pcl_from_plyfile(std::string plyfile_location){
    happly::PLYData plyIn(plyfile_location);
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

    return pclmatrix;
}

int main(int argc, char* argv[]) {

    using namespace std;
    // load the ply file
    if(argc != 3){
        cout << "Kindly enter the location of the plyfile of the modelview and dataview\n";
        return 1;
    }

    char* modelview_location = argv[1];
    char* dataview_location = argv[2];

    Eigen::Matrix dataview = get_pcl_from_plyfile(dataview_location);
    Eigen::Matrix modelview = get_pcl_from_plyfile(modelview_location);

    // calculate centroid of model view
    Eigen::Matrix<float, 3, 1> model_mean = calculate_centroid(modelview);
    //cout << model_mean;
    // calculate centroid of data view
    Eigen::Matrix<float, 3, 1> data_mean = calculate_centroid(dataview);

    float sse = 1.1; // arbitrary value initialization
    int iteration = 0;

    // get the time now
    auto start = chrono::steady_clock::now();
    while (sse > 0.1) {
        // get correspondent points
        Eigen::MatrixXf correspondent_points(3, dataview.cols());
        correspondent_points = get_correspondence_points(modelview, dataview);
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

        int determinant_prod = det_u * det_v;
        if (determinant_prod == 1) {
            s.diagonal() << 1, 1, 1;
        } else if (determinant_prod == -1) {
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
        Eigen::Matrix<float, 4, 4> transformational_matrix = get_transformation_matrix(rotational_matrix, translational_matrix);
        // cout << "\n transformational_matrix is " << transformational_matrix << "\n";
        dataview = transform_matrix(dataview, rotational_matrix, translational_matrix);
        data_mean = calculate_centroid(dataview);
        sse = calculate_sse(dataview, modelview);

        cout << "\nSSE is " << sse << " in iteration " << iteration++;
    }
    auto end = chrono::steady_clock::now();
    auto time_taken = chrono::duration_cast<chrono::milliseconds>(end-start).count();

    cout << "\nTotal time taken for the ICP process is " << time_taken << " milliseconds";

}
