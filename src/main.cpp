// Prevent the from printing the default instructions for use
#ifndef IGL_VIEWER_VIEWER_QUIET
#define IGL_VIEWER_VIEWER_QUIET
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <vector>
#include <chrono>
#include <igl/readOFF.h>
#include <igl/arap.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/unproject.h>
#include <igl/project.h>
#include "rotation_matrix.h"
#include "cotangent_weight_matrix.h"
#include "system_matrix.h"
#include "ARAP_iteration.h"
#include "main.h"

void updateMesh(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd &U, const Eigen::MatrixXi &F, const std::set<int> &fixedPoints, const Eigen::RowVector3d &pointColor, int currentHandle, const Eigen::RowVector3d &handleColor)
{
    viewer.data().clear();
    viewer.data().clear_points();
    viewer.data().set_mesh(U, F);

    // draw the fixed points
    for (int pointIndex : fixedPoints)
    {
        viewer.data().add_points(U.row(pointIndex), pointColor);
    }

    // draw the handle point
    viewer.data().add_points(U.row(currentHandle), handleColor);
    viewer.core().align_camera_center(U, F);
    viewer.data().set_face_based(true);
}

void extractScreenCoords(igl::opengl::glfw::Viewer viewer, double &x, double &y)
{
    x = viewer.current_mouse_x;
    y = viewer.core().viewport(3) - viewer.current_mouse_y;
}

Eigen::Vector3d getWorldCoords(igl::opengl::glfw::Viewer viewer)
{
    double x, y;
    extractScreenCoords(viewer, x, y);

    Eigen::Vector3d worldProjection = igl::unproject(Eigen::Vector3f(x, y, 0), viewer.core().view,
                                                     viewer.core().proj, viewer.core().viewport)
                                          .cast<double>();

    return worldProjection;
}

Eigen::Vector3d getWorldCoords(igl::opengl::glfw::Viewer viewer, Eigen::Vector3d handle)
{
    double x, y;
    extractScreenCoords(viewer, x, y);
    Eigen::Vector3f handleFloat = handle.cast<float>();
    Eigen::Vector3f meshProjection = igl::project(handleFloat, viewer.core().view, viewer.core().proj, viewer.core().viewport);

    Eigen::Vector3d worldProjection = igl::unproject(Eigen::Vector3f(x, y, meshProjection(2)), viewer.core().view,
                                                     viewer.core().proj, viewer.core().viewport)
                                          .cast<double>();

    // std::cout << "world " << worldProjection.x() << " " << worldProjection.y() << " " << worldProjection.z() << std::endl;

    return worldProjection;
}

void saveFixedPoints(const std::set<int> &fixedPoints, const std::string &mesh)
{
    std::fstream file("../fixedPoints.txt", std::fstream::out | std::fstream::trunc);

    std::cout << "Saving points" << std::endl;
    for (int point : fixedPoints)
    {
        file << point << std::endl;
        // std::cout << point << std::endl;
    }

    file.close();
}

void loadFixedPoints(std::set<int> &fixedPoints)
{
    std::fstream file("../fixedPoints.txt", std::fstream::in);

    std::string point;
    std::cout << "Loading points" << std::endl;

    while (std::getline(file, point))
    {
        fixedPoints.insert(std::stoi(point, nullptr));
        // std::cout << point << std::endl;
    }
    file.close();
}

void ourImplArap(Eigen::Vector3d projection)
{
    bool FIRST_LOOP = true;
    for (int iteration = 0; iteration < MAX_ITER; iteration++)
    {
        if (FIRST_LOOP)
        {
            U = V;

            U.row(currentHandle) = projection;

            FIRST_LOOP = false;
        }
        // 1st:
        //  Calculate the vector R of rotation matrices
        //  U is the OUT matrix of the previous iteration
        std::vector<Eigen::Matrix3d> R = rotation_matrix(PD, U, neighbors);

        // Note: we need to make an initial guess for U at the first step! See paper to understand why?
        // 2nd:
        // U is the matrix of the new vertices after one ARAP iteration, called OUT or U above
        U = ARAP_iteration(fixedPoints, W, R, V, neighbors, L_initial);
        std::cout << "Done with arap: " << iteration << std::endl;
        // Note: Inside the ARAP_iteration function the matrix L_initial is copied to a matrix called L since it will be modified
        // (constraints addition or removal), but the L_initial will stay the same during the process
        // Repeat 1 and 2 until satisfaction(Energy difference < epsilon or number of iterations or both)
    }

    FIRST_LOOP = true;
}

void libiglImplArap(Eigen::Vector3d projection, igl::ARAPEnergyType energy = igl::ARAP_ENERGY_TYPE_SPOKES)
{
    // Use the ARAP implementation from the igl library
    igl::ARAPData arap_data;
    arap_data.max_iter = MAX_ITER;
    arap_data.with_dynamics = true;
    arap_data.energy = energy;
    // arap_data.energy = igl::ARAP_ENERGY_TYPE_SPOKES;
    // arap_data.energy = igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS;
    // arap_data.energy = igl::ARAP_ENERGY_TYPE_ELEMENTS;

    // Convert the set<int> to a vector a Eigen::VectorXi of positions
    Eigen::VectorXi fixedPointsEigen = Eigen::Map<Eigen::VectorXi>(std::vector<int>(fixedPoints.begin(), fixedPoints.end()).data(), fixedPoints.size());

    // Add the handle
    fixedPointsEigen.conservativeResize(fixedPointsEigen.size() + 1);
    fixedPointsEigen(fixedPointsEigen.size() - 1) = currentHandle;

    // Create the matrix of fixed points, (.size() is also counting the handle now)
    Eigen::MatrixXd FixedPointsMatrix(fixedPointsEigen.size(), V.cols());
    // Populate the matrix
    for (unsigned int i = 0; i < fixedPointsEigen.size() - 1; ++i)
    {
        FixedPointsMatrix.row(i) = V.row(fixedPointsEigen[i]);
    }
    // Add the handle
    FixedPointsMatrix.row(FixedPointsMatrix.rows() - 1) = V.row(currentHandle);

    FixedPointsMatrix.row(FixedPointsMatrix.rows() - 1) = projection;

    // Precomputation
    igl::arap_precomputation(V, F, V.cols(), fixedPointsEigen, arap_data);

    // Solve
    igl::arap_solve(FixedPointsMatrix, arap_data, U);
}

double calculateEnergy()
{
    double energy = 0;
    int vertices = V.rows();
    std::vector<Eigen::Matrix3d> R = rotation_matrix(PD, U, neighbors);
    for (int i = 0; i < vertices; i++)
    {
        for (int j = 0; j < neighbors[i].size(); j++)
        {
            energy += W(i, j) * ((U.row(i) - U.row(j)).transpose() - R[i] * (V.row(i) - V.row(j)).transpose()).norm();
        }
    }

    return energy;
}

void testLibiglArap(Eigen::Vector3d proj, igl::ARAPEnergyType energy)
{
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "MAX iteration " << MAX_ITER << " with energy " << energy << std::endl;
    libiglImplArap(proj, energy);
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "It took " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000000.0 << " seconds." << std::endl;
    std::cout << "The calculated energy is " << calculateEnergy() << std::endl;
}

void testOurArap(Eigen::Vector3d proj)
{
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "MAX iteration " << MAX_ITER << std::endl;
    ourImplArap(proj);
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "It took " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000000.0 << " seconds." << std::endl;
    std::cout << "The calculated energy is " << calculateEnergy() << std::endl;
}

void runTests()
{
    currentHandle = 850;
    Eigen::Vector3d proj(0.792619, 0.0459771, -0.222237);
    useIglArap = true;
    std::cout << "Running with libigl's implementation:\n";
    igl::readOFF(mesh, V, F);
    loadFixedPoints(fixedPoints);
    MAX_ITER = 1;
    testLibiglArap(proj, igl::ARAP_ENERGY_TYPE_SPOKES);
    testLibiglArap(proj, igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS);
    testLibiglArap(proj, igl::ARAP_ENERGY_TYPE_ELEMENTS);
    MAX_ITER = 5;
    testLibiglArap(proj, igl::ARAP_ENERGY_TYPE_SPOKES);
    testLibiglArap(proj, igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS);
    testLibiglArap(proj, igl::ARAP_ENERGY_TYPE_ELEMENTS);
    MAX_ITER = 10;
    testLibiglArap(proj, igl::ARAP_ENERGY_TYPE_SPOKES);
    testLibiglArap(proj, igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS);
    testLibiglArap(proj, igl::ARAP_ENERGY_TYPE_ELEMENTS);
    MAX_ITER = 20;
    testLibiglArap(proj, igl::ARAP_ENERGY_TYPE_SPOKES);
    testLibiglArap(proj, igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS);
    testLibiglArap(proj, igl::ARAP_ENERGY_TYPE_ELEMENTS);

    std::cout << "\n\n\nRunning with our implementation:\n";
    MAX_ITER = 1;
    testOurArap(proj);
    MAX_ITER = 5;
    testOurArap(proj);
    MAX_ITER = 10;
    testOurArap(proj);
    MAX_ITER = 20;
    testOurArap(proj);
}

int main()
{

    // Print our instructions to the console
    std::cout << Instructions << std::endl;

    while (true)
    {
        std::cout << "Please choose a mesh:\n1 - Stanford Bunny\n2 - Armadillo\n";
        std::cin >> meshChoice;
        if (meshChoice == "1")
        {
            mesh = bunnyMesh;
            break;
        }
        else if (meshChoice == "2")
        {
            mesh = armadilloMesh;
            break;
        }
        else if (meshChoice == "3")
        {
            mesh = exMesh;
            break;
        }
        else
        {
            std::cout << "The choice is invalid.\n";
        }
    }
    igl::readOFF(mesh, V, F);

    unsigned int vertices = V.rows();
    neighbors = std::vector<std::vector<unsigned int>>(vertices);
    // Compute the neighbors list
    igl::adjacency_list(F, neighbors);
    // Initialize matrix U
    U = Eigen::MatrixXd(V.rows(), V.cols());
    // Compute the weight matrix for the mesh
    W = weight_matrix(V, F, neighbors);
    // Compute the constant part
    PD = compute_const_part_covariance(V, W, neighbors);
    // Save the original system matrix
    L_initial = initialize_system_matrix(V, W, neighbors);

    while (true)
    {
        std::cout << "Please choose an ARAP implementation:\n1 - Our implementation\n2 - igl library implementation\n3 - Run tests\n";
        std::cin >> arapChoice;
        if (arapChoice == "1")
        {
            useIglArap = false;
            break;
        }
        else if (arapChoice == "2")
        {
            useIglArap = true;
            break;
        }
        else if (arapChoice == "3")
        {
            runTests();
            return 0;
        }
        else
        {
            std::cout << "The choice is invalid.\n";
        }
    }

    viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer &viewer, int button, int) -> bool
    {
        if (mode != ROTATION)
        {
            // Some of the code in this function is taken from the official libigl tutorial at https://github.com/libigl/libigl/blob/main/tutorial/708_Picking/main.cpp
            // the face id of the face that was clicked on
            int fid;
            // the barycentric coordinates of the face being clicked on
            Eigen::Vector3f bc;

            double x, y;
            extractScreenCoords(viewer, x, y);

            if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
                                         viewer.core().proj, viewer.core().viewport, V, F, fid, bc))
            {

                Eigen::Vector3i selectedFace = F.row(fid);
                // maxCoeff will give the point index whose barycentric coordinate is the greatest meaning that is the closest one to the actual point that was clicked
                // Eigen::Vector3f point = V(selectedFace(bc.maxCoeff()));
                int closestPointIndex;
                bc.maxCoeff(&closestPointIndex);
                int pointIndex = selectedFace(closestPointIndex);

                bool remove = (button == RIGHT_BUTTON);

                if (mode == HANDLE_SELECTION)
                {
                    currentHandle = pointIndex;
                    viewer.data().add_points(V.row(pointIndex), handleColor);
                }
                else if (mode == POINT_SELECTION)
                {
                    if (remove)
                    {
                        fixedPoints.erase(pointIndex);
                        viewer.data().clear_points();
                        for (auto &point : fixedPoints)
                        {
                            viewer.data().add_points(V.row(point), pointColor);
                        }
                    }
                    else
                    {
                        std::cout << "Inserting a fixed point: " << pointIndex << "  At coordinate: " << V.row(pointIndex) << std::endl;
                        fixedPoints.insert(pointIndex);
                        viewer.data().add_points(V.row(pointIndex), pointColor);
                    }
                }

                return true;
            }
        }

        return false;
    };

    viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer &viewer, int button, int) -> bool
    {
        bool remove = (button == RIGHT_BUTTON);

        if (mode == HANDLE_SELECTION && currentHandle != -1)
        {
            Eigen::Vector3d projection = getWorldCoords(viewer, V.row(currentHandle));
            if (useIglArap)
            {
                libiglImplArap(projection);
            }
            else
            {
                ourImplArap(projection);
            }

            updateMesh(viewer, U, F, fixedPoints, pointColor, currentHandle, handleColor);
            currentHandle = -1;
            return true;
        }
        return false;
    };

    viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int) -> bool
    {
        std::cout << "Pressed key " << key << std::endl;
        switch (key)
        {
        case 'H':
        case 'h':
            mode = HANDLE_SELECTION;
            std::cout << "Switching to handle" << std::endl;
            return true;
        case 'F':
        case 'f':
            mode = POINT_SELECTION;
            std::cout << "Switching to point" << std::endl;
            return true;
        case 'N':
        case 'n':
            mode = ROTATION;
            std::cout << "Switching to rotation" << std::endl;
            return true;
        case 'R':
        case 'r':
            igl::readOFF(mesh, V, F);
            fixedPoints.clear();
            viewer.data().clear();
            viewer.data().clear_points();
            viewer.data().set_face_based(true);
            viewer.data().set_mesh(V, F);
            return true;
        case 's':
        case 'S':
            saveFixedPoints(fixedPoints, mesh);
            return true;
        case 'l':
        case 'L':
            loadFixedPoints(fixedPoints);
            updateMesh(viewer, V, F, fixedPoints, pointColor, currentHandle, handleColor);
            return true;
        }
        return true;
    };

    viewer.data().set_mesh(V, F);
    viewer.core().align_camera_center(V, F);
    viewer.data().set_face_based(true);
    viewer.data().point_size = 10;
    viewer.launch();
}