#pragma once

#include <Eigen/Core>
#include <vector>
#include <igl/adjacency_list.h>

 /**
  * Performs a single iteration of the ARAP algorithm.
  * @param constraint_vertices matrix containing the vertices that make up the constraints of the problem
  * @param constraint_indices vector of indices for the vertices that make up the contraints of the problem
  * @param W matrix of weights for the entire mesh (1s in the main diagonal)  
  * @param R vector of rotation matricies for every vertex i
  * @param V matrix where each row represents a vertex coordinates in 3D space
  * @param L_initial matrix of weights "system matrix" from paper
  * @return OUT matrix where each row represents a vertex coordinates in 3D space after it has been moved
  */
Eigen::MatrixXd ARAP_iteration(
    const Eigen::MatrixXd& constraint_vertices,
    const Eigen::VectorXi& constraint_indices,
    const Eigen::MatrixXd& W,
    const std::vector<Eigen::Matrix3d>& R,
    const Eigen::MatrixXd& V,
    const std::vector<std::vector<unsigned int>>& neighbors,
    const Eigen::MatrixXd& L_initial);