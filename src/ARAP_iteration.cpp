#include "ARAP_iteration.h"

Eigen::MatrixXd ARAP_iteration(
  const std::set<int>& constraint_indices,
  const Eigen::MatrixXd& W,
  const std::vector<Eigen::Matrix3d>& R,
  const Eigen::MatrixXd& V,	
  const std::vector<std::vector<unsigned int>>& neighbors,
  const Eigen::MatrixXd& L_initial)
{
    // Amount of vertices
	unsigned int vertices = V.rows();		

    // Make a copy so L can be modified by the constraints
    Eigen::MatrixXd L = L_initial;

    // Initialize the right hand side of the system
    Eigen::VectorXd b_x = Eigen::VectorXd::Zero(vertices);
    Eigen::VectorXd b_y = Eigen::VectorXd::Zero(vertices);
    Eigen::VectorXd b_z = Eigen::VectorXd::Zero(vertices);

    // Compute the right hand side 
    for (unsigned int i = 0; i < vertices; ++i) {
        // sum vector
        Eigen::Vector3d sum (0.0, 0.0, 0.0);         
        // loop over the one-ring neighborhood of the current vertex i
		// for each vertex, calculate the sum over its one-ring neighborhood, considering cotangent weights and rotation matrices
        for (unsigned int j: neighbors[i]) {
			// for each neighboring vertex j of vertex i, calculate a weight w_ij from the cotangent weight matrix W
			// w_ij is associated with the edge between the vertices i and j
            double w_ij = W.coeff(i, j); 
			// compute the local contribution to the sum (see eq. 8 from the reference article)
            //std::cout << i << " " << j << " " << vertices;
            //std::cout << "second" <<R.size() << " " << V.rows() << " " << vertices;
            sum += w_ij/2.0*(R[i]+R[j])*(V.row(i)-V.row(j)).transpose(); 
        }
		// the local contributions for all neighboring vertices are accumulated in the sum vector and assigned to the corresponding 
        // components of b for the current vertex i
        b_x(i) = sum.x();
        b_y(i) = sum.y();
        b_z(i) = sum.z();
    }

    // Adjust the b vectors based on user-defined constraints (see eq. 10 from the reference article)
	// for each constraint, modify the entries in the b vectors, taking into account the system matrix entries corresponding to the constrained vertices
	// update the system matrix by zeroing out rows and columns corresponding to the constrained vertices and set the diagonal entry to 1.0

    // loop over constraints
    for (const auto& c : constraint_indices) {
	
        // for each constraint, retrieve the constraint vector and the corresponding index from the respective matrices
        Eigen::Vector3d constraint = V.row(c);
        int constraint_index = c;

        for (unsigned int i = 0; i < vertices; ++i) {
            //for each vertex i check if the system matrix coefficient at position (i, constraint_index) is not zero
            if (L.coeff(i, constraint_index) != 0.0) {      

			// adjust the b vectors by subtracting the product of the system matrix coefficient and the corresponding constraint component
			// the (i, constraint_index) position corresponds to the interaction between the vertex i and the constrained degree of freedom represented by constraint_index
                b_x(i) -= L.coeff(i, constraint_index)*constraint.x();
                b_y(i) -= L.coeff(i, constraint_index)*constraint.y();
                b_z(i) -= L.coeff(i, constraint_index)*constraint.z();
            }
        }
		// set the value of the constraint itself in the b vectors at the constraint's index 
        b_x(constraint_index) = constraint.x();
        b_y(constraint_index) = constraint.y();
        b_z(constraint_index) = constraint.z();

        // modify the system matrix to account for the constraints
		// ensure that the constrained degrees of freedom are fixed at a specifiec value, and the corresponding entries in the b vectors are set to the constrained values
		// this modification effectively freezes the constrained degree of freedom, preventing it from contributing to the deformation
        for (unsigned int i = 0; i < vertices; ++i) {
            if(L.coeff(constraint_index, i) != 0.0) { 
                //set the coefficient at position (constraint_index, i) and (i, constraint_index) to zero
                L(constraint_index, i) = 0.0; 
                L(i, constraint_index) = 0.0;
                L(constraint_index, constraint_index) = 1.0; 								 
            }
        }
    }

    // use Cholesky solver for new vertex positions
    Eigen::LLT<Eigen::MatrixXd> solver(L);
	// cholesky factorization is a numerical method for solving symmetric positive definite linear systems, and it decomposes 
    // the system matrix into the product of a lower triangular matrix and its transpose
    Eigen::VectorXd x = solver.solve(b_x);
    Eigen::VectorXd y = solver.solve(b_y);
    Eigen::VectorXd z = solver.solve(b_z);

    // create the output matrix for the now moved vertices
    Eigen::MatrixXd OUT (vertices, 3);

    // pupulate it  
    for (unsigned int i = 0; i < vertices; ++i) {
        OUT.row(i) = Eigen::Vector3d(x(i), y(i), z(i)); 
    }

    return OUT;
}
