//
// Created by brazenparadox on 29/08/23.
//

#include<eigen3/Eigen/Dense>

Eigen::MatrixXf rotate_matrix_along_z(Eigen::MatrixXf matrix, float degree){
    float theta = degree * 3.14/180;
    Eigen::Matrix<float, 3, 3> zrotational_matrix;
    zrotational_matrix << cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1 ;
    return zrotational_matrix*matrix;
}

Eigen::MatrixXf rotate_matrix_along_x(Eigen::MatrixXf matrix, float degree){
    float theta = degree * 3.14/180;
    Eigen::Matrix<float, 3, 3> xrotational_matrix;
    xrotational_matrix << 1, 0, 0, 0, cos(theta), -sin(theta), 0, sin(theta), cos(theta);
    return xrotational_matrix*matrix;
}

Eigen::MatrixXf rotate_matrix_along_y(Eigen::MatrixXf matrix, float degree){
    float theta = degree * 3.14/180;
    Eigen::Matrix<float, 3, 3> yrotational_matrix;
    yrotational_matrix << cos(theta), 0, -sin(theta),  0, 1, 0, sin(theta), 0, cos(theta)  ;
    return yrotational_matrix*matrix;
}

Eigen::MatrixXf translate_matrix(Eigen::MatrixXf matrix, float x, float y, float z){
    Eigen::Vector3f translation_matrix = Eigen::Vector3f(x, y, z);
    //cout << translation_matrix;
    return matrix.colwise() + translation_matrix;
}

Eigen::Matrix<float, 4, 4> get_transformation_matrix(Eigen::Matrix<float, 3, 3> rotational_matrix, Eigen::Matrix<float, 3, 1> translational_matrix){

    Eigen::Matrix<float, 4, 4> transformational_matrix;
    for(int row=0; row<3; row++){
        for(int col=0; col<3; col++){
            transformational_matrix(row, col) = rotational_matrix(row, col);
        }
    }

    for(int row=0; row<3; row++){
        transformational_matrix(row, 3) = translational_matrix(row, 0);
    }

    for(int col=0; col<3; col++){
        transformational_matrix(3, col) = 0;
    }

    transformational_matrix(3, 3) = 1;

    return transformational_matrix;
}
