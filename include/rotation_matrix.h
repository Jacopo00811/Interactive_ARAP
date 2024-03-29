#pragma once

#include <vector>
#include <Eigen/Core>
#include <igl/adjacency_list.h>


/* 
The rotation matrix is found by SVD of the covariance matrix S. The covariance matrix S is made up by a constant part P*D
and a non constant part P'^(T). Where P and P'^(T) is a 3 x |N(v_i)| containing e_ij's or e_ij''s as its colunmns and D is the
cotangent weight matrix.
*/

/**
 * Computes the constant part (P*D) of the covariance matrix S given:
 * @param V matrix where each row represents a vertex coordinates in 3D space
 * @param W matrix of weights for the entire mesh (1s in the main diagonal)
 * @return vector of matricies of PD
 */
std::vector<Eigen::Matrix<double, 3, -1>> compute_const_part_covariance(const Eigen::MatrixXd& V, const Eigen::MatrixXd& W, const std::vector<std::vector<unsigned int>>& neighbors);

/**
 * Computes a vector of rotation matricies R given:
 * @param PD constant part of the covariance matrix S
 * @param OUT matrix where each row represents a vertex coordinates in 3D space after it has been moved
 * @return vector of matricies of R
 */
std::vector<Eigen::Matrix3d> rotation_matrix(const std::vector<Eigen::Matrix<double, 3, -1>>& PD, const Eigen::MatrixXd& OUT, const std::vector<std::vector<unsigned int>>& neighbors);
