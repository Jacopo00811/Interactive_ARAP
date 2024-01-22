#include "one_iteration.h"
#include "cotangent_weight_matrix.h"
#include "rotation_matrix.h"
#include <igl/min_quad_with_fixed.h> //igl library used for solving quadratic optimization problems (minimizing energy) with fixed variables (mesh subject to fixed value boundary conditions)
#include <igl/adjacency_list.h>      //igl library which constructs the graph adjacency list of a given mesh (v,f)
									 //returns a list containing at index 'i' the adjacent vertices (neighbors) of vertex 'i'

void one_iteration(
  const Eigen::MatrixXd& constr_vertices,
  const Eigen::VectorXi& constr_indices,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& U, 							
  Eigen::SparseMatrix<double>& systemMatrix_init, 
  bool unif_weight_use,
  double unif_weight_val)
{
	unsigned int numVertex = U.rows(); //no. of vertices
	unsigned int numFaces = F.rows();  //no. of faces

	//**********Local Optimization

	//std::vector<Eigen::Matrix3d> R (numVertex); 					//rotation matrix
	//std::vector<Eigen::Matrix<double, 3, -1>> P_prime(numVertex); // -1 indicates that the no. of columns can vary at runtime


    std::vector<std::vector<double>> neighbors(numVertex); //defines a vector of neighbors for one vertex
    igl::adjacency_list(F, neighbors); 					   //populate the 'neighbors' vector with the neighbors of each vertex 'i'

/*// Compose 'P_prime' with the 'e_prime_ij'
    for (unsigned int i = 0; i < numVertex; ++i) { 
        for (unsigned int j : neighbors[i]) { 
            Eigen::Vector3d e_prime_ij = U.row(i) - U.row(j); 						 //difference between the positions of a vertex and its neighbors
            P_prime[i].conservativeResize(P_prime[i].rows(), P_prime[i].cols() + 1); //resize the matrix to rows x cols while leaving old values untouched (for matrices of dynamic size)
            P_prime[i].col(P_prime[i].cols() - 1) = e_prime_ij; 					 //add 'e_prime_ij' to 'P_prime' as the last column
        }
    }*/

	/*//solve for each rotation matrix
	for (int i = 0; i < numVertex; i++) {
		//calculate the covariance matrix for each vertex
		Eigen::Matrix3d S_i = K[i] * P_prime[i].transpose();
//perform Singular Value Decomposition (SVD) on the covariance matrix using the Eigen::JacobiSVD function
//SVD is a matrix factorization method that decomposes a matrix into three separate matrices: U, Σ, and V. It is used to solve linear equations and find the rotation matrices in this code.
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(S_i,Eigen::ComputeFullU|Eigen::ComputeFullV);
		//calculate the rotation matrix using the SVD results
		R[i] = svd.matrixV() * svd.matrixU().transpose();
	}*/

	//**********Global Optimization
	//*****System Initialization
    Eigen::SparseMatrix<double> systemMatrix = systemMatrix_init; //copy 'systemMatrix_init' to 'systemMatrix'
	//define the Right-Hand Side vectors (see eq. 8 in the reference paper)
    Eigen::VectorXd m_rhsX(numVertex);
    Eigen::VectorXd m_rhsY(numVertex);
    Eigen::VectorXd m_rhsZ(numVertex);

    //*****Initialize the Right-Hand Side (RHS) vectors 
    m_rhsX.setZero();
    m_rhsY.setZero();
    m_rhsZ.setZero();

    //*****Calculate the cotangent matrix W (with uniform weights) 'w_ij'
    Eigen::SparseMatrix<double> W;
    compute_cotangent_weight_for_pair(V, F, W, unif_weight_use,unif_weight_val);

    //*****Compute the RHS
    // Loop over each vertex
    for (unsigned int i = 0; i < numVertex; ++i) {
        Eigen::Vector3d sum(0.0, 0.0, 0.0);         //define and initialize the sum vector

        //Loop over the one-ring neighborhood of the current vertex 'i'
		//For each vertex, calculate the sum over its one-ring neighborhood, considering cotangent weights and rotation matrices
        for (unsigned int j: neighbors[i]) {

			//for each neighboring vertex 'j' of vertex 'i', calculate a weight 'w_ij' from the cotangent matrix 'W'
			//'w_ij' is associated with the edge between the vertices 'i' and 'j'
            double w_ij = W.coeff(i, j); 

			//local contribution to the sum using the corangent weight 'w_ij', rotation matrices 'R[i]' and 'R[j]' and the difference in vertex positions ('V.row(i)-V.row(j)')
			//see eq 8 from the reference article		
            sum += w_ij / 2.0 * (R[i] + R[j]) * (V.row(i) - V.row(j)).transpose(); 
        }
		//the computed local contributions for all neighboring vertices are accumulated in the 'sum' vector and assigned to the corresponding components of the RHS vectors for the current vertex 'i'
        m_rhsX(i) = sum.x();
        m_rhsY(i) = sum.y();
        m_rhsZ(i) = sum.z();
    }

    //*****Incorporate constraints (see eq. 10 from the reference article)
	//Adjust the RHS vectors based on user-defined constraints
	//For each constraint, modify the entries in the RHS vectors, taking into account the system matrix entries corresponding to the constrained vertices
	//Update the system matrix by zeroing out rows and columns corresponding to the constrained vertices and set the diagonal entry to 1.0

    for (int constraint = 0; constraint < constr_vertices.rows(); constraint++) { //loop over constraints
	//for each constraint, retrieve the constraint vector 'c' and the corresponding index 'c_index' from the respective matrices
        Eigen::Vector3d c = constr_vertices.row(constraint);
        int c_index = constr_indices(constraint);

        // Adapt RHS
        for (unsigned int i = 0; i < numVertex; ++i) {
            if (systemMatrix.coeff(i, c_index) != 0.0) {      //for each vertex 'i' check if the system matrix/s coefficient at position '(i,c_index)' is not zero

			//adjust the RHS vectors (m_rhsX, m_rhsY, m_rhsZ) by subtracting the product of the system matrix coefficient and the corresponding constraint component.
			//the '(i, c_index)' position corresponds to the interaction between the vertex 'i' and the constrained degree of freedom represented by 'c_index'.
                m_rhsX(i) -= systemMatrix.coeff(i, c_index) * c.x();
                m_rhsY(i) -= systemMatrix.coeff(i, c_index) * c.y();
                m_rhsZ(i) -= systemMatrix.coeff(i, c_index) * c.z();
            }
        }
		//set the value of the constraint itself in the RHS vectors at the constraint's index 'c_index'.
        m_rhsX(c_index) = c.x();
        m_rhsY(c_index) = c.y();
        m_rhsZ(c_index) = c.z();

        //*****Delete rows and columns
		//modify the system matrix to account for the constraints
		//ensure that the constrained degrees of freedom are fixed, and the corresponding entries in the RHS vectors are set to the constrained values
		//This modification effectively freezes the constrained degree of freedom, preventing it from contributing to the deformation.

        for (unsigned int i = 0; i < numVertex ; ++i) { 
            if(systemMatrix.coeff(c_index, i) != 0.0) { 
                systemMatrix.coeffRef(c_index, i) = systemMatrix.coeffRef(i, c_index) = 0.0; //set the coefficient at position '(c_index, i)' and '(i, c_index)' to zero.
                systemMatrix.coeffRef(c_index, c_index) = 1.0; 								 //if 'i' is equal to 'c_index', set the coefficient at position '(c_index, c_index)' to 1.0.
				//the diagonal entry at '(c_index, c_index)' corresponds to the constrained degree of freedom itself within the system matrix
				//ensures that the constrained degree of freedom remains fixed at its specified value
            }
        }
    }

    //*****Solve for new vertex positions
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> svd(systemMatrix);
	//Eigen::SimplicialCholesky is a solver class in the Eigen library that performs the Cholesky factorization on the sparse system matrix ('systemMatrix').
	//Cholesky factorization is a numerical method for solving symmetric positive definite linear systems, and it decomposes the system matrix into the product of a lower triangular matrix and its transpose.
    Eigen::VectorXd x = svd.solve(m_rhsX);
    Eigen::VectorXd y = svd.solve(m_rhsY);
    Eigen::VectorXd z = svd.solve(m_rhsZ);

    for (unsigned int i = 0; i < numVertex; ++i) {
        U.row(i) = Eigen::Vector3d(x(i), y(i), z(i));    //for each vertex 'i', set the new position in the matrix 'U' (which contains all deformed vertices).
		//The x, y, and z components of the new position are taken from the solutions obtained in the linear system solving step (x(i), y(i), z(i)).
    }

}
