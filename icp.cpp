//
// Created by brazenparadox on 25/8/23.
//

#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/Jacobi>
#include "iostream"


Eigen::Matrix<float, 3, 1> calculate_centroid(Eigen::MatrixXf pcl){
    Eigen::Matrix<float, 3, 1> centroid;
    long column_size = pcl.cols();

    for(int row=0; row<3; row++){
        float sum = 0;
        for(int col=0; col<column_size; col++){
            sum += pcl(row, col);
        }
        centroid(row, 0) = sum/column_size;
    }

    return centroid;
}

Eigen::MatrixXf  get_correspondence_points(Eigen::MatrixXf modelview, Eigen::MatrixXf dataview){
    /*This methods calculates the euclidean distance between the points in modelview and dataview pointset.
     * A new matrix is returned with the model points that are corresponding to the data points in their respective indices */

    Eigen::MatrixXf correspondent_points(3, dataview.cols());
    for(int column=0; column<dataview.cols(); column++){
        float smallest_distance;
        for(int model_column=0; model_column<modelview.cols(); model_column++){
            float euc_dist = sqrt(pow(modelview(0, model_column) - dataview(0, column), 2) +
                                           pow(modelview(1, model_column) - dataview(1, column), 2) +
                                           pow(modelview(2, model_column) - dataview(2, column), 2));
            if (model_column == 0 || euc_dist < smallest_distance){
                smallest_distance = euc_dist;
                for(int i=0; i<3; i++)
                    correspondent_points(i, column) = modelview(i, model_column);
            }
        }
    }
    return correspondent_points;
}

Eigen::MatrixXf get_meaned_data(Eigen::MatrixXf data, Eigen::Matrix<float, 3, 1> mean){
    for(int col=0; col<data.cols(); col++){
        for(int row=0; row<3; row++)
            data(row, col) = data(row, col) - mean(row, 0);
    }
    return data;
}

Eigen::Matrix<float, 3, 3> calculate_cross_covariance(Eigen::MatrixXf dataview, Eigen::MatrixXf correspondences){
    Eigen::Matrix<float, 3, 1> data_mean = calculate_centroid(dataview);
    Eigen::Matrix<float, 3, 1> correspondence_mean = calculate_centroid(correspondences);

    Eigen::MatrixXf meaned_correspondence = get_meaned_data(correspondences, correspondence_mean);
    //std::cout << "\n\nThe meaned correspondent points are " << meaned_correspondence;
    Eigen::MatrixXf meaned_data = get_meaned_data(dataview, data_mean);

//    std::cout << "\n  correspondences : " << correspondences;
//    std::cout << "\n meaned correspondence : " << meaned_correspondence << "\n\n meaned data: " << meaned_data;

    Eigen::Matrix<float, 3, 3> cross_covariance_matrix = meaned_data*meaned_correspondence.transpose();
//    std::cout << cross_covariance_matrix.rows() << "<-- rows   columns --> " << cross_covariance_matrix.cols();

    // Round the values of the covariance matrix
    for(int row=0; row<cross_covariance_matrix.rows(); row++){
        for(int col=0; col<cross_covariance_matrix.cols(); col++){
            float temp = round(cross_covariance_matrix(row, col) * 1000);
            cross_covariance_matrix(row, col) = temp/1000;
        }
    }

    return cross_covariance_matrix;
}

Eigen::MatrixXf transform_matrix(Eigen::MatrixXf matrix, Eigen::Matrix<float, 3, 3> rotational_matrix, Eigen::Matrix<float, 3, 1> translation_matrix){
    Eigen::MatrixXf transformed_matrix = rotational_matrix * matrix;
    transformed_matrix = transformed_matrix.colwise() + translation_matrix;
    return transformed_matrix;
}

float calculate_sse(Eigen::MatrixXf dataview, Eigen::MatrixXf correspondent_points){
    float sse = 0;
    float euc_dist = 0;
    for(int column=0; column<dataview.cols(); column++){
         euc_dist = sqrt(pow(correspondent_points(0, column) - dataview(0, column), 2) +
                              pow(correspondent_points(1, column) - dataview(1, column), 2) +
                              pow(correspondent_points(2, column) - dataview(2, column), 2));
        sse += euc_dist;
    }
    return sse;
}
