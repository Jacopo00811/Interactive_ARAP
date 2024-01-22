#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include "cotangent_weight_matrix.h"
#include "rotation_matrix.h"
#include "one_iteration.h"

int main(int argc, char* argv[])
{
    const Eigen::MatrixXd V{
    {1.0, 0.0, 2.0},
    {2.0, 0.0, 2.0},
    {2.0, 2.0, 2.0},
    {2.0, -3.0, 2.0},
    {4.0, 0.0, 2.0}
    };

    Eigen::MatrixXd U{
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

    auto PD = compute_const_part_covariance(V, F);
    std::cout << std::endl << std::endl << std::endl;
    std::cout << "PD (weights * e_ij): " << std::endl;

    for (const auto& matrix : PD) {
    	std::cout << "Matrix: " << std::endl << matrix << std::endl << std::endl;
    }

    auto R = rotation_matrix(PD, F, U);
    std::cout << std::endl << std::endl << std::endl;
    std::cout << "Rotation matricies: " << std::endl;

    for (const auto& matrix : R) {
        std::cout << "Matrix: " << std::endl << matrix << std::endl << std::endl;
    }
}