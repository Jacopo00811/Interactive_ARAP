#include <system_matrix.h>

Eigen::MatrixXd initialize_system_matrix(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& W) {
	// Amount of vertices
	unsigned int vertices = V.rows();

	Eigen::MatrixXd L(vertices, vertices);

	std::vector<std::vector<unsigned int>> neighbors(vertices);
	igl::adjacency_list(F, neighbors);

	for (unsigned int i = 0; i < vertices; ++i) {
		for (unsigned int j : neighbors[i]) {
			L(i, i) += W.coeff(i, j);
			L(i, j) = -W.coeff(i, j);
		}
	}
	return L;
}
