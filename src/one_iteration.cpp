#include "one_iteration.h"
#include "cotangent_weight_matrix.h"
#include "rotation_matrix.h"
#include <igl/min_quad_with_fixed.h> //igl library used for solving quadratic optimization problems (minimizing energy) with fixed variables (mesh subject to fixed value boundary conditions)
#include <igl/adjacency_list.h>      //igl library which constructs the graph adjacency list of a given mesh (v,f)
									 //returns a list containing at index 'i' the adjacent vertices (neighbors) of vertex 'i'

void one_iteration(
  const Eigen::MatrixXd& constr_vertices,
  const Eigen::VectorXi& constr_indices,
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& P_prime, 							
  Eigen::SparseMatrix<double>& systemMatrix_init)
{
	unsigned int numVertex = P_prime.rows();               //no. of vertices
	unsigned int numFaces  = F.rows();                     //no. of faces
    std::vector<std::vector<double>> neighbors(numVertex); //defines a vector of neighbors for one vertex
    igl::adjacency_list(F, neighbors); 					   //populate the 'neighbors' vector with the neighbors of each vertex 'i'

	//*****Initialize the system
    //use a local variable to work with the initial system matrix
    Eigen::SparseMatrix<double> systemMatrix = systemMatrix_init; 
	//define the Right-Hand Side vectors (see eq. 8 in the reference paper)
    Eigen::VectorXd rhs_x(numVertex);
    Eigen::VectorXd rhs_y(numVertex);
    Eigen::VectorXd rhs_z(numVertex);
    //initialize the Right-Hand Side (RHS) vectors with zero
    rhs_x.setZero();
    rhs_y.setZero();
    rhs_z.setZero();

    //*****Compute the RHS vectors
    for (unsigned int i = 0; i < numVertex; ++i) {
        //define and initialize the sum vector
        Eigen::Vector3d sum(0.0, 0.0, 0.0);         
        //loop over the one-ring neighborhood of the current vertex 'i'
		//for each vertex, calculate the sum over its one-ring neighborhood, considering cotangent weights (see 'cotangent_weight_matrix.cpp') and rotation matrices (see 'rotation_matrix.cpp')
        for (unsigned int j: neighbors[i]) {
			//for each neighboring vertex 'j' of vertex 'i', calculate a weight 'w_ij' from the cotangent weight matrix 'W'
			//'w_ij' is associated with the edge between the vertices 'i' and 'j'
            double w_ij = W.coeff(i, j); 
			//compute the local contribution to the sum using the cotangent weight 'w_ij', rotation matrices 'R[i]' and 'R[j]' and the difference in vertex positions ('P.row(i)-P.row(j)')
			//see eq 8 from the reference article		
            sum += w_ij / 2.0 * (R[i] + R[j]) * (P.row(i) - P.row(j)).transpose(); 
        }
		//the computed local contributions for all neighboring vertices are accumulated in the 'sum' vector and assigned to the corresponding components of the RHS vectors for the current vertex 'i'
        rhs_x(i) = sum.x();
        rhs_y(i) = sum.y();
        rhs_z(i) = sum.z();
    }

    //*****Adjust the RHS vectors based on user-defined constraints (see eq. 10 from the reference article)
	//For each constraint, modify the entries in the RHS vectors, taking into account the system matrix entries corresponding to the constrained vertices
	//Update the system matrix by zeroing out rows and columns corresponding to the constrained vertices and set the diagonal entry to 1.0

    //loop over constraints
    for (int c = 0; c < constr_vertices.rows(); c++) {
	//for each constraint, retrieve the constraint vector 'constraint' and the corresponding index 'constraint_in' from the respective matrices
        Eigen::Vector3d constraint = constr_vertices.row(c);
        int constraint_in = constr_indices(c);
        for (unsigned int i = 0; i < numVertex; ++i) {
            //for each vertex 'i' check if the system matrix coefficient at position '(i,constraint_in)' is not zero
            if (systemMatrix.coeff(i, constraint_in) != 0.0) {      
			//adjust the RHS vectors (rhs_x, rhs_y, rhs_z) by subtracting the product of the system matrix coefficient and the corresponding constraint component.
			//the '(i, constraint_in)' position corresponds to the interaction between the vertex 'i' and the constrained degree of freedom represented by 'constraint_in'.
                rhs_x(i) -= systemMatrix.coeff(i, constraint_in) * constraint.x();
                rhs_y(i) -= systemMatrix.coeff(i, constraint_in) * constraint.y();
                rhs_z(i) -= systemMatrix.coeff(i, constraint_in) * constraint.z();
            }
        }
		//set the value of the constraint itself in the RHS vectors at the constraint's index 'constraint_in'.
        rhs_x(constraint_in) = constraint.x();
        rhs_y(constraint_in) = constraint.y();
        rhs_z(constraint_in) = constraint.z();

        //*****Modify the system matrix to account for the constraints
		//Ensure that the constrained degrees of freedom are fixed at a specifiec value, and the corresponding entries in the RHS vectors are set to the constrained values
		//This modification effectively freezes the constrained degree of freedom, preventing it from contributing to the deformation.

        for (unsigned int i = 0; i < numVertex ; ++i) { 
            if(systemMatrix.coeff(constraint_in, i) != 0.0) { 
                systemMatrix.coeffRef(constraint_in, i) = systemMatrix.coeffRef(i, constraint_in) = 0.0; //set the coefficient at position '(constraint_in, i)' and '(i, constraint_in)' to zero.
                systemMatrix.coeffRef(constraint_in, constraint_in) = 1.0; 								 //if 'i' is equal to 'constraint_in', set the coefficient at position '(constraint_in, constraint_in)' to 1.0.
				//the diagonal entry at '(constraint_in, constraint_in)' corresponds to the constrained degree of freedom itself within the system matrix
            }
        }
    }

    //*****Use Cholesky solver for new vertex positions
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> svd(systemMatrix);
	//Eigen::SimplicialCholesky is a solver class in the Eigen library that performs the Cholesky factorization on the sparse system matrix ('systemMatrix').
	//Cholesky factorization is a numerical method for solving symmetric positive definite linear systems, and it decomposes the system matrix into the product of a lower triangular matrix and its transpose.
    Eigen::VectorXd x = svd.solve(rhs_x);
    Eigen::VectorXd y = svd.solve(rhs_y);
    Eigen::VectorXd z = svd.solve(rhs_z);

    for (unsigned int i = 0; i < numVertex; ++i) {
        P_prime.row(i) = Eigen::Vector3d(x(i), y(i), z(i));    
        //for each vertex 'i', set the new position in the matrix 'P_prime' (which contains all deformed vertices).
		//the x, y, and z components of the new position are taken from the solutions obtained in the linear system solving step (x(i), y(i), z(i)).
    }

}
