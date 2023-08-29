//
// Created by brazenparadox on 25/8/23.
//

#ifndef INC_3D_RESIGTRATION_CASTELLANI_AND_BARTOLI_CPP_ICP_H
#define INC_3D_RESIGTRATION_CASTELLANI_AND_BARTOLI_CPP_ICP_H

#endif //INC_3D_RESIGTRATION_CASTELLANI_AND_BARTOLI_CPP_ICP_H

#include<eigen3/Eigen/Dense>


Eigen::Matrix<float, 3, 1> calculate_centroid(Eigen::MatrixXf);
Eigen::MatrixXf get_correspondence_points(Eigen::MatrixXf, Eigen::MatrixXf);
Eigen::Matrix<float, 3, 3> calculate_cross_covariance(Eigen::MatrixXf, Eigen::MatrixXf);
Eigen::MatrixXf transform_matrix(Eigen::MatrixXf, Eigen::Matrix<float, 3, 3>, Eigen::Matrix<float, 3, 1>);
float calculate_sse(Eigen::MatrixXf, Eigen::MatrixXf);
//int calculate_rotation(v, s, u);
//int calculate_translation(model_mean, rotational_matrix, data_mean);
