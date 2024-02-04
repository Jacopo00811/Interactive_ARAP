#pragma once

void updateMesh(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& U, const Eigen::MatrixXi& F, const std::set<int>& fixedPoints, const Eigen::RowVector3d& pointColor, int currentHandle, const Eigen::RowVector3d& handleColor);
