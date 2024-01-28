#include <Eigen/Core>   //includes the main module of Eigen provding dense matrix and vector support
#include <Eigen/Sparse> //includes sparse matrices in Eigen 
                        //a sparse matrix is a matrix in which majority of the elements are 0
/**
 * 'one_iteration.h' defines a void function named 'one_iteration' that performs a single iteration of the ARAP algorithm. 
 * The function takes several input parameters and updates the 'U' matrix and 'systemMatrix_init' based on the deformation.
 * For a const reference, the function body cannot directly change the value of that object
 * 
 * param K                    = A vector of matrices representing the bi-linear term of energy combining rotations and positions.
 *                              The bi-linear term of energy combines rotations and positions to measure the rigidity of the mesh. 
 *                              It is pre-constructed.
 * 
 * @param constr_vertices     = A matrix containing the control point positions. (constraint vertices)
 *                                  Control points are specific vertices in the mesh that are constrained to certain positions. 
 * @param constr_indices      = A vector containing the indices of the control points in the vertex matrix 'V' (indices of constraint vertices)
 * 
 * @param P                   = Matrix containing all vertices in the initial position
 * @param F                   = Matrix containin all the faces of the mesh
 * @param P_prime             = Matrix containing all deformed vertices (the positions of the vertices after deformation). 
 *                              This matrix will be updated by the function.
 * 
 * @param systemMatrix_init   = A sparse matrix representing the system matrix for global optimization. 
 *                              It represents the linear system of equations to be solved during the deformation process.
 *                              This matrix will be updated by the function.
 * 
 */

void one_iteration(
  const std::vector<Eigen::Matrix<double,3,-1>>& K,
  const Eigen::MatrixXd& constr_vertices,
  const Eigen::VectorXi& constr_indices,P,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& P_prime,
  Eigen::SparseMatrix<double>& systemMatrix_init);