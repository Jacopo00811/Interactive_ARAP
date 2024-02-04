#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <igl/readOFF.h>
#include <igl/arap.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/unproject.h>
#include "rotation_matrix.h"
#include "cotangent_weight_matrix.h"
#include "system_matrix.h"
#include "ARAP_iteration.h"
#include "main.h"

enum Mode {
	ROTATION,
	POINT_SELECTION,
	HANDLE_SELECTION
};

// similar to the code in https://libigl.github.io/dox/Viewer_8h_source.html
enum MouseButton {
	LEFT_BUTTON,
	MIDDLE_BUTTON,
	RIGHT_BUTTON
};

void updateMesh(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& U, const Eigen::MatrixXi& F, const std::set<int>& fixedPoints, const Eigen::RowVector3d& pointColor, int currentHandle, const Eigen::RowVector3d& handleColor)
{
    viewer.data().clear();
    viewer.data().clear_points();
    viewer.data().set_mesh(U, F);

    // draw the fixed points
    for (int pointIndex : fixedPoints) {
        viewer.data().add_points(U.row(pointIndex), pointColor);
    }

    // draw the handle point
    viewer.data().add_points(U.row(currentHandle), handleColor);
    viewer.core().align_camera_center(U, F);
    viewer.data().set_face_based(true);
}


int main()
{
    bool useIglArap;

    const std::string pathToMeshes = "../Meshes/";
    const std::string bunnyMesh = pathToMeshes + "bunny.off";
    const std::string armadilloMesh = pathToMeshes + "armadillo.off";
    const std::string exMesh = pathToMeshes + "ex.off";

    std::string meshChoice;
    std::string mesh;

    bool FIRST_LOOP = true;
    int MAX_ITER = 5;


    Mode mode = Mode::ROTATION;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::RowVector3d handleColor(0, 255, 0);
    Eigen::RowVector3d pointColor(0, 0, 255);
  
    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;

    int currentHandle {0};
    std::set<int> fixedPoints; // better to be kept as indices rather than concrete points
    std::string arapChoice;

    // Find the neighbors of each vertex
    std::vector<Eigen::Matrix<double, 3, -1>> PD;
    Eigen::MatrixXd W;
    Eigen::MatrixXd L_initial;

    const std::string Instructions = R"(
    [left click]                        Place fixed point when in fixed point mode
    [left click] + [drag]               Place and move handle point when in handle point mode
    [left click] (not on mesh)          Rotation
    H,h                                 Enter handle point mode
    F,f                                 Enter fixed point mode
    N,n                                 Enter rotation mode (default at the start)
    )";
    
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

    while (true)
    {
        std::cout << "Please choose an ARAP implementation:\n1 - Our implementation\n2 - igl library implementation\n";
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
        else
        {
            std::cout << "The choice is invalid.\n";
        }
    }

    unsigned int vertices = V.rows();
    std::vector<std::vector<unsigned int>> neighbors(vertices);
    // Compute the neighbors list
    igl::adjacency_list(F, neighbors);
    // Initialize matrix U
    Eigen::MatrixXd U(V.rows(), V.cols());
    // Compute the weight matrix for the mesh
    W = weight_matrix(V, F, neighbors);
    // Compute the constant part 
    PD = compute_const_part_covariance(V, W, neighbors);
    // Save the original system matrix
    L_initial = initialize_system_matrix(V, W, neighbors);

    viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer &viewer, int button, int) -> bool
    {
        // Some of the code in this function is taken from the official libigl tutorial at https://github.com/libigl/libigl/blob/main/tutorial/708_Picking/main.cpp
        // the face id of the face that was clicked on
        int fid;
        // the barycentric coordinates of the face being clicked on
        Eigen::Vector3f bc;

        double x = viewer.current_mouse_x;
        double y = viewer.core().viewport(3) - viewer.current_mouse_y;
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

        return false;
    };

    viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer &viewer, int button, int) -> bool
    {
        bool remove = (button == RIGHT_BUTTON);

        if (mode == HANDLE_SELECTION)
        {
            if (useIglArap)
            {
                // Use the ARAP implementation from the igl library
                igl::ARAPData arap_data;
                arap_data.max_iter = MAX_ITER;
                arap_data.with_dynamics = true;

                // Convert the set<int> to a vector a Eigen::VectorXi of positions
                Eigen::VectorXi fixedPointsEigen = Eigen::Map<Eigen::VectorXi>(std::vector<int>(fixedPoints.begin(), fixedPoints.end()).data(), fixedPoints.size());
                
                // Add the handle
                fixedPointsEigen.conservativeResize(fixedPointsEigen.size()+1);
                fixedPointsEigen(fixedPointsEigen.size()-1) = currentHandle;

                // Create the matrix of fixed points, (.size() is also counting the handle now)
                Eigen::MatrixXd FixedPointsMatrix(fixedPointsEigen.size(), V.cols());
                // Populate the matrix
                for (unsigned int i = 0; i < fixedPointsEigen.size()-1; ++i) {
                    FixedPointsMatrix.row(i) = V.row(fixedPointsEigen[i]);
                }
                // Add the handle
                FixedPointsMatrix.row(FixedPointsMatrix.rows()-1) = V.row(currentHandle);

                // Update the position of the current handle based on the mouse position
                double x = viewer.current_mouse_x;
                double y = viewer.core().viewport(3) - viewer.current_mouse_y;
                Eigen::Vector3d projection = igl::unproject(Eigen::Vector3f(x, y, 0), viewer.core().view,
                    viewer.core().proj, viewer.core().viewport).cast<double>();
                FixedPointsMatrix.row(FixedPointsMatrix.rows()-1) -= projection;

                // Precomputation
                igl::arap_precomputation(V, F, V.cols(), fixedPointsEigen, arap_data);

                // Solve
                igl::arap_solve(FixedPointsMatrix, arap_data, U);
                
                updateMesh(viewer, U, F, fixedPoints, pointColor, currentHandle, handleColor);

                return true;
            }
            else
            {
                for (int iteration = 0; iteration < MAX_ITER; iteration++) {
                    if (FIRST_LOOP) {
                        U = V;

                        double x = viewer.current_mouse_x;
                        double y = viewer.core().viewport(3) - viewer.current_mouse_y;
                        Eigen::Vector3d projection = igl::unproject(Eigen::Vector3f(x, y, 0), viewer.core().view,
                            viewer.core().proj, viewer.core().viewport).cast<double>();
                  
                        U.row(currentHandle) -= projection;

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
                updateMesh(viewer, U, F, fixedPoints, pointColor, currentHandle, handleColor);

                return true;
            }  
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
        }
        return true;
    };

    viewer.data().set_mesh(V, F);
    viewer.core().align_camera_center(V, F);
    viewer.data().set_face_based(true);
    viewer.launch();
}