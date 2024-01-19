#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/adjacency_list.h>

/**
 * Given two vectors of coordinates, computes the cotangent of the angle in between them.
 * @param a coordinates of first vector
 * @param b coordinates of second vector
 * @return cot(angle) between vector a and b
 */
double cotangent(Eigen::Vector3d a, Eigen::Vector3d b);

/**
 * Given two indices of two vetecies in V and a face (i.e. vector of indices corresponding to a row of matrix F), returns true if
 * both the indicies are in the face (triangle).
 * @param i first index
 * @param j first index
 * @param face face (vector of int)
 * @return true if i and j are in f
 */
bool is_part_of_face(unsigned int i, unsigned int j, Eigen::Vector3i face);

/**
 * Given two indices of vetecies in V and a face (i.e. vector of indices corresponding to a row of matrix F), returns the index of 
 * the third vertex.
 * @param i first index
 * @param j first index
 * @param face face (vector of int)
 * @return index of third vertex
 */
int get_third_face_vertex(unsigned int i, unsigned int j, Eigen::Vector3i face);

/**
 * Constructs the cotagent matrix D of weights given:
 * @param V matrix where each row represents a vertex coordinates in 3D space 
 * @param F matrix where each row represents the three vertices that compose a face f (triangle)
 * @return W matrix of weights for the entire mesh (1s in the main diagonal)  
 */
Eigen::MatrixXd weight_matrix(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

/**
 * Helper function, computes the cotanget weight for a pair of indices
 * @param i first index
 * @param j first index
 * @param V matrix where each row represents a vertex coordinates in 3D space
 * @param F matrix where each row represents the three vertices that compose a face f (triangle)
 * @return weight_ij 
 */
double compute_cotangent_weight_for_pair(unsigned int i, unsigned int j, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);