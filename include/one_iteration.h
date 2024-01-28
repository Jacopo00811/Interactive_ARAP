#include <Eigen/Core>   //includes the main module of Eigen provding dense matrix and vector support
#include <Eigen/Sparse> //includes sparse matrices in Eigen 
                        //a sparse matrix is a matrix in which majority of the elements are 0
/**
 * 'one_iteration.h' defines a void function named 'one_iteration' that performs a single iteration of the ARAP algorithm. 
 * The function takes several input parameters and returns the 'U' matrix based on the deformation.
 * For a const reference, the function body cannot directly change the value of that object
 * 
 * @param constr_vertices a matrix containing the control point positions. (constraint vertices)
 *                        Control points are specific vertices in the mesh that are constrained to certain positions. 
 * @param constr_indices a vector containing the indices of the control points in the vertex matrix 'V' (indices of constraint vertices)
 * @param W weight matrix computed in 'cotangent_weight_matrix.cpp'
 * @param R rotation matrix computed in 'rotation_matrix.cpp'
 * 
 * @param P matrix containing all vertices in the initial position
 * @param F matrix containin all the faces of the mesh
 * @param P_prime matrix containing all deformed vertices (the positions of the vertices after deformation). Corresponds to 'U' in 'rotation_matrix.cpp'
 *                              This matrix will be updated and returned by the function.
 * 
 * @param systemMatrix_init a sparse matrix representing the system matrix for global optimization. 
 *                              It represents the linear system of equations to be solved during the deformation process.
 *                              This matrix will be updated by the function.
 * 
 */

Eigen::MatrixXd one_iteration(
  const Eigen::MatrixXd& constr_vertices,
  const Eigen::VectorXi& constr_indices,
  const Eigen::MatrixXd& W,
  const std::vector<Eigen::Matrix3d>& R,
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& P_prime,
  Eigen::SparseMatrix<double>& systemMatrix_init);