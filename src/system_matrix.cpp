#include "system_matrix.h"

Eigen::MatrixXd initialize_system_matrix(const Eigen::MatrixXd& V, const Eigen::MatrixXd& W, const std::vector<std::vector<unsigned int>>& neighbors) {
	// Amount of vertices
	unsigned int vertices = V.rows();

	Eigen::MatrixXd L(vertices, vertices);

	for (unsigned int i = 0; i < vertices; ++i) {
		for (unsigned int j : neighbors[i]) {
			L(i, i) += W.coeff(i, j);
			L(i, j) = -W.coeff(i, j);
		}
	}
	return L;
}
