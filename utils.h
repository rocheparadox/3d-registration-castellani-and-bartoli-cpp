//
// Created by brazenparadox on 29/08/23.
//

#ifndef INC_3D_RESIGTRATION_CASTELLANI_AND_BARTOLI_CPP_UTILS_H
#define INC_3D_RESIGTRATION_CASTELLANI_AND_BARTOLI_CPP_UTILS_H

#endif //INC_3D_RESIGTRATION_CASTELLANI_AND_BARTOLI_CPP_UTILS_H

#include<eigen3/Eigen/Dense>

Eigen::MatrixXf rotate_matrix_along_z(Eigen::MatrixXf, float);
Eigen::MatrixXf rotate_matrix_along_x(Eigen::MatrixXf, float);
Eigen::MatrixXf rotate_matrix_along_y(Eigen::MatrixXf, float);
Eigen::MatrixXf translate_matrix(Eigen::MatrixXf, float, float, float);
Eigen::Matrix<float, 4, 4> get_transformation_matrix(Eigen::Matrix<float, 3, 3>, Eigen::Matrix<float, 3, 1>);