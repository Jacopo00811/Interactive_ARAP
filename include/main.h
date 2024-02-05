#pragma once

void updateMesh(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& U, const Eigen::MatrixXi& F, const std::set<int>& fixedPoints, const Eigen::RowVector3d& pointColor, int currentHandle, const Eigen::RowVector3d& handleColor);

/**
 * Calculates the screen coordinates based on the position of the mouse
 * @param viewer the viewer object
 * @param x the result for the x coordinate
 * @param y the result for the y coordinate
*/
void extractScreenCoords(igl::opengl::glfw::Viewer viewer, double& x, double& y);

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