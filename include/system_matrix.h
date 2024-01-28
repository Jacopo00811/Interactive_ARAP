#pragma once
#include <Eigen/Core>
#include <igl/adjacency_list.h>

/**
 * Initializes the system matrix L so it can be used in in the ARAP iteration.
 * @param V matrix where each row represents a vertex coordinates in 3D space
 * @param W matrix of weights for the entire mesh (1s in the main diagonal)
 * @param neighbors vector of vectors of neighbors
 * @return system matrix L
 */
Eigen::MatrixXd initialize_system_matrix(const Eigen::MatrixXd& V, const Eigen::MatrixXd& W, const std::vector<std::vector<unsigned int>>& neighbors);