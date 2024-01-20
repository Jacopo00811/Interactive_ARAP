#include "cotangent_weight_matrix.h"

double cotangent(Eigen::Vector3d a, Eigen::Vector3d b) {
	return (a.dot(b)) / (a.cross(b)).norm();
}

bool is_part_of_face(unsigned int i, unsigned int j, Eigen::Vector3i face) {
	return (face[0] == i && (face[1] == j || face[2] == j)) ||
		(face[1] == i && (face[0] == j || face[2] == j)) ||
		(face[2] == i && (face[1] == j || face[0] == j));
}

int get_third_face_vertex(unsigned int i, unsigned int j, Eigen::Vector3i face) {
	for (unsigned int idx : face) {
		if (idx != i && idx != j) {
			return idx;
		}
	}
	return -1;
}

Eigen::MatrixXd weight_matrix(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
	// Amount of vertices
	unsigned int vertices = V.rows();
	
	// Create the weight matrix
	Eigen::MatrixXd W = Eigen::MatrixXd::Zero(vertices, vertices);
	
	// Find the neighbors of each vertex
	std::vector<std::vector<unsigned int>> neighbors(vertices);
	igl::adjacency_list(F, neighbors);
	// For vertex i, neighbors[i] contains the indices of all the vertices that are connected to vertex i
	for (unsigned int i = 0; i < vertices; i++) {
		
		// Find all row indices (faces) that contain the target value i (vertex i)
			// Step 1: Create a binary matrix indicating where the target value i is present
			// and find all row indices (faces) that contain the target value i (vertex i)
		Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> mask = (F.array() == i);

		std::vector<Eigen::Vector3i> faces_with_vertex_i;
    
		for (int k = 0; k < mask.rows(); ++k) {
			// Check if the row has at least one occurrence of the target value i
			if (mask.row(k).any()) {
				// Add the row
				faces_with_vertex_i.push_back(F.row(k));
			}
		}

		for (unsigned int j : neighbors[i]) {
			
			double w_ij = 0;
			
			if (W.coeff(j, i) == 0)
				w_ij = compute_cotangent_weight_for_pair(i, j, faces_with_vertex_i, V, F);
			else
				w_ij = W.coeff(j, i);

			W(i, j) = w_ij;
		}
		W(i, i) = 1;
	}

	return W;
}

double compute_cotangent_weight_for_pair(unsigned int i, unsigned int j, const std::vector<Eigen::Vector3i>& faces_with_vertex_i, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
	// To store the faces that contain the vertex i 
	std::vector<Eigen::Vector3i> possible_faces;

	// Find all possible faces that contain both vertex i and vertex j
	for (Eigen::Vector3i face : faces_with_vertex_i) { 
		if (is_part_of_face(i, j, face))
			possible_faces.push_back(face);
	}
	// Check the size of the vector
	assert(possible_faces.size() <= 2);
	// Get values of vertices i and j
	Eigen::Vector3d v_i = V.row(i);
	Eigen::Vector3d v_j = V.row(j);

	// Compute: 0.5*(cot(alpha)+cot(betha))
	double cotangent_sum = 0.0f;
	// Compute the cotangent for each face and save result in cotangent_sum
	for (Eigen::Vector3i face : possible_faces) {
		int vertex_id = get_third_face_vertex(i, j, face);
		Eigen::Vector3d v_o = V.row(vertex_id);
		double angle = cotangent(v_i - v_o, v_j - v_o);
		cotangent_sum += angle;
	}
	// Divide by two (1/2*abs((cot(alpha)+cot(betha)))
	return abs(cotangent_sum)/2.0f;
}