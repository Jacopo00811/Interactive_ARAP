#pragma once

enum Mode
{
    ROTATION,
    POINT_SELECTION,
    HANDLE_SELECTION
};

// similar to the code in https://libigl.github.io/dox/Viewer_8h_source.html
enum MouseButton
{
    LEFT_BUTTON,
    MIDDLE_BUTTON,
    RIGHT_BUTTON
};

void updateMesh(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd &U, const Eigen::MatrixXi &F, const std::set<int> &fixedPoints, const Eigen::RowVector3d &pointColor, int currentHandle, const Eigen::RowVector3d &handleColor);

/**
 * Calculates the screen coordinates based on the position of the mouse
 * @param viewer the viewer object
 * @param x the result for the x coordinate
 * @param y the result for the y coordinate
 */
void extractScreenCoords(igl::opengl::glfw::Viewer viewer, double &x, double &y);

/**
 * Calculates the projection of the coordinates based position of the mouse in the world space
 * @param viewer the viewer object
 * @return the projection vector
 */
Eigen::Vector3d getWorldCoords(igl::opengl::glfw::Viewer viewer);

/**
 * Calculates the projection of the coordinates based position of the mouse and the handle in the world space
 * @param viewer the viewer object
 * @param handle the coordinates of the handle
 * @return the projection vector
 */
Eigen::Vector3d getWorldCoords(igl::opengl::glfw::Viewer viewer, Eigen::Vector3d handle);

bool useIglArap;

const std::string pathToMeshes = "../Meshes/";
const std::string bunnyMesh = pathToMeshes + "bunny.off";
const std::string armadilloMesh = pathToMeshes + "armadillo.off";
const std::string exMesh = pathToMeshes + "ex.off";

std::string meshChoice;
std::string mesh;

bool FIRST_LOOP = true;
int MAX_ITER = 1;

Mode mode = Mode::ROTATION;

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::RowVector3d handleColor(0, 255, 0);
Eigen::RowVector3d pointColor(0, 0, 255);

// Plot the mesh
igl::opengl::glfw::Viewer viewer;

int currentHandle{0};
std::set<int> fixedPoints; // better to be kept as indices rather than concrete points
std::string arapChoice;

// Find the neighbors of each vertex
std::vector<Eigen::Matrix<double, 3, -1>> PD;
Eigen::MatrixXd W;
Eigen::MatrixXd L_initial;
Eigen::MatrixXd U;
std::vector<std::vector<unsigned int>> neighbors;

const std::string Instructions = R"(
[left click]                        Place fixed point when in fixed point mode
[left click] + [drag]               Place and move handle point when in handle mode , or rotate when in rotation mode
[right click]                       Delete a point
H,h                                 Enter handle point mode
F,f                                 Enter fixed point mode
N,n                                 Enter rotation mode (default at the start)
R,r                                 Reset the mesh
)";