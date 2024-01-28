#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include "cotangent_weight_matrix.h"
#include "rotation_matrix.h"
#include "system_matrix.h"
#include "ARAP_iteration.h"

int main(int argc, char* argv[])
{
    const Eigen::MatrixXd V{
    {1.0, 0.0, 2.0},
    {2.0, 0.0, 2.0},
    {2.0, 2.0, 2.0},
    {2.0, -3.0, 2.0},
    {4.0, 0.0, 2.0}
    };

    const Eigen::MatrixXd U{
    {2.0, 0.0, 2.0},
    {10.0, 0.0, 2.0},
    {9.0, 6.0, 8.0},
    {7.0, -3.0, 3.0},
    {2.0, 0.0, 1.0}
    };

    const Eigen::MatrixXi F{
            {0, 1, 2},
            {0, 1, 3},
            {1, 2, 4},
            {1, 3, 4},
    };

    const Eigen::MatrixXd constraint_vertices{
    {6.4, 2.1, 3.1},
    {1.0, 2.0, 3.0},
    };

    // Create a vector of indices in this weird way    
    const Eigen::VectorXi constraint_indices = [] {
        Eigen::VectorXi indices(3);
        indices << 1, 3, 4;
        indices << 0, 1, 2;
        return indices;
        }();

    /*
    // Inline mesh of a cube
    const Eigen::MatrixXd V = (Eigen::MatrixXd(8, 3) <<
        0.0, 0.0, 0.0,
        0.0, 0.0, 1.0,
        0.0, 1.0, 0.0,
        0.0, 1.0, 1.0,
        1.0, 0.0, 0.0,
        1.0, 0.0, 1.0,
        1.0, 1.0, 0.0,
        1.0, 1.0, 1.0).finished();
    const Eigen::MatrixXi F = (Eigen::MatrixXi(12, 3) <<
        0, 6, 4,
        0, 2, 6,
        0, 3, 2,
        0, 1, 3,
        2, 7, 6,
        2, 3, 7,
        4, 6, 7,
        4, 7, 5,
        0, 4, 5,
        0, 5, 1,
        1, 5, 7,
        1, 7, 3).finished();
    
    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
    viewer.launch();
    */



    // Amount of vertices
    unsigned int vertices = V.rows();
    // Find the neighbors of each vertex
    std::vector<std::vector<unsigned int>> neighbors(vertices);
    igl::adjacency_list(F, neighbors);
    auto W = weight_matrix(V, F, neighbors);

    auto PD = compute_const_part_covariance(V, W, neighbors);
    std::cout << std::endl << std::endl << std::endl;
    std::cout << "PD (weights * e_ij): " << std::endl;
    for (const auto& matrix : PD) {
    	std::cout << "Matrix: " << std::endl << matrix << std::endl << std::endl;
    }

    auto R = rotation_matrix(PD, U, neighbors);
    std::cout << std::endl << std::endl << std::endl;
    std::cout << "Rotation matricies: " << std::endl;
    for (const auto& matrix : R) {
        std::cout << "Matrix: " << std::endl << matrix << std::endl << std::endl;
    }
    
    auto L = initialize_system_matrix(V, W, neighbors);
    std::cout << std::endl << std:: endl << "Matrix L: " << std::endl << L << std::endl << std::endl;
    
    auto res = ARAP_iteration(constraint_vertices, constraint_indices, W, R, V, neighbors, L);
    std::cout << "Matrix res: " << std::endl << res << std::endl << std::endl;
    std::cout << "Matrix V: " << std::endl << V << std::endl << std::endl;
    std::cout << std::endl << std::endl << "Matrix L after: " << std::endl << L << std::endl << std::endl;

    /*
    Steps to perform one ARAP itaration:
        
        Initialize:
            // Amount of vertices
            unsigned int vertices = V.rows();
            
            // Find the neighbors of each vertex
            std::vector<std::vector<unsigned int>> neighbors(vertices);
            igl::adjacency_list(F, neighbors);
            
            auto W = weight_matrix(V, F, neighbors);
        
            auto PD = compute_const_part_covariance(V, W, neighbors);

            // Save the original system matrix
            auto L_initial = initialize_system_matrix(V, W, neighbors);
        

        1st:
            // Calculate the vector R of rotation matrices
            // U is the OUT matrix of the previous iteration
            auto R = rotation_matrix(PD, U, neighbors);
            
            // Note: we need to make an initial guess for U at the first step! See paper to understand why?
        2nd:
            // res is the matrix of the new vertices after one ARAP iteration, called OUT or U above
            auto res = ARAP_iteration(constraint_vertices, constraint_indices, W, R, V, neighbors, L_initial);

            // Note: Inside the ARAP_iteration function the matrix L_initial is copied to a matrix called L since it will be modified 
            // (constraints addition or removal), but the L_initial will stay the same during the process

        Repeat 1 and 2 until satisfaction (Energy difference < epsilo or number of iterations or both)
    */

}