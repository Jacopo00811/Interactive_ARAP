#pragma once
#include <Eigen/Core>
#include <igl/adjacency_list.h>

/**
 * Initializes the system matrix L so it can be used in in the ARAP iteration.
 * @param V matrix where each row represents a vertex coordinates in 3D space
 * @param F matrix where each row represents the three vertices that compose a face f (triangle)
 * @param F W matrix of weights for the entire mesh (1s in the main diagonal)  
 * @return system matrix L
 */
Eigen::MatrixXd initialize_system_matrix(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& W);