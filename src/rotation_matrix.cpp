#include <rotation_matrix.h>
#include <cotangent_weight_matrix.h>

std::vector<Eigen::Matrix<double, 3, -1>> compute_const_part_covariance(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
	// Amount of vertices
	unsigned int vertices = V.rows();

	std::vector<Eigen::Matrix<double, 3, -1>> PD(vertices);
	//std::vector<Eigen::MatrixXd> D(vertices);

	Eigen::MatrixXd W = weight_matrix(V, F);

	std::vector<std::vector<unsigned int>> neighbors(vertices);
	igl::adjacency_list(F, neighbors);

	// For vertex i, neighbors[i] and the weights compute the matrix PD[i]
	for (unsigned int i = 0; i < vertices; ++i) {
		//Eigen::MatrixXd D_i = Eigen::MatrixXd::Zero(vertices, vertices);
		for (unsigned int j : neighbors[i]) {	
			Eigen::Vector3d e_ij = V.row(i) - V.row(j);
			PD[i].conservativeResize(PD[i].rows(), PD[i].cols() + 1);
			//D_i(j, j) = W.coeff(i, j);
			PD[i].col(PD[i].cols()-1) = W.coeff(i, j) * e_ij;            
		}
	}


	return PD;
}

std::vector<Eigen::Matrix3d> rotation_matrix(const std::vector<Eigen::Matrix<double, 3, -1>>& PD, const Eigen::MatrixXi& F, Eigen::MatrixXd& U) {
	// Amount of vertices
	unsigned int vertices = U.rows();

	// Vectors of rotation matricies R and P_prime
	std::vector<Eigen::Matrix3d> R(vertices);
	std::vector<Eigen::Matrix<double, 3, -1>> P_prime(vertices);

	std::vector<std::vector<double>> neighbors(vertices);
	igl::adjacency_list(F, neighbors);
	// Compose P_prime with the e_prime_ij
	for (unsigned int i = 0; i < vertices; ++i) {
		for (unsigned int j : neighbors[i]) {
			Eigen::Vector3d e_prime_ij = U.row(i) - U.row(j);
			P_prime[i].conservativeResize(P_prime[i].rows(), P_prime[i].cols() + 1);
			P_prime[i].col(P_prime[i].cols() - 1) = e_prime_ij;
		}
	}


	// Solve for each rotation matrix R_i in the vector R
	for (int i = 0; i < vertices; i++) {
		// Calculate each covariance matrix S_i
		Eigen::Matrix3d S_i = PD[i] * P_prime[i].transpose();
		// SVD decomposition on S_i to find V and U^(T)
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(S_i, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::Matrix3d temp = svd.matrixV() * svd.matrixU().transpose();
		// Check and correct determinant less than zero
		if (temp.determinant() < 0) {
			Eigen::MatrixXd SVD_U = svd.matrixU();
			SVD_U.rightCols(1) = SVD_U.rightCols(1) * -1;
			temp = svd.matrixV() * SVD_U.transpose();
		}
		assert(temp.determinant() > 0);
		// Update the vector of rotation matricies
		R[i] = temp;
	}

	return R;
}