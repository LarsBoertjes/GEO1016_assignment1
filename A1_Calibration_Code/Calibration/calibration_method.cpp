/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "calibration.h"
#include "matrix_algo.h"


using namespace easy3d;

std::vector<Vector3D> points_3d = {{5, 0, 2}, {3, 8, 0}, {0, 7, 7}, {1, 0, 1},
                                   {10, 10, 0}, {5, 3, 0}, {1, 9, 0}, {3, 0, 8},
                                   {10, 0, 10}, {2, 8, 0}};

std::vector<Vector2D> points_2d = {{494, 464}, {473, 239}, {297, 301}, {414, 435},
                                   {662, 204}, {518, 374}, {426, 207}, {365, 524},
                                   {536, 610}, {450, 237}};

/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 * @return True on success, otherwise false. On success, the camera parameters are returned by fx, fy, cx, cy, skew, R, and t).
 */
bool Calibration::calibration(
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx,  /// output: focal length (i.e., K[0][0]).
        double& fy,  /// output: focal length (i.e., K[1][1]).
        double& cx,  /// output: x component of the principal point (i.e., K[0][2]).
        double& cy,  /// output: y component of the principal point (i.e., K[1][2]).
        double& s,   /// output: skew factor (i.e., K[0][1]), which is s = -alpha * cot(theta).
        Matrix33& R, /// output: the 3x3 rotation matrix encoding camera rotation.
        Vector3D& t) /// output：a 3D vector encoding camera translation.
{

    // TODO: check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    if (points_3d.size() != points_2d.size() || points_3d.size() < 6) {
        std::cerr << "Error: The number of 3D and 2D points must be equal and more than 6" << std::endl;
        return false;
    }

    // TODO: construct the P matrix (so P * m = 0).
    // Done using the lecture notes 02-camera_calibration paragraph 2
    int num_points = points_3d.size();
    Matrix P(2 * num_points, 12);

    for (int i = 0; i < num_points; ++i) {
        const Vector3D& X = points_3d[i];  // X = 3D point
        const Vector2D& x = points_2d[i];  // x = corresponding 2D point

        // X-coordinate row
        P(2 * i, 0) = X.x();
        P(2 * i, 1) = X.y();
        P(2 * i, 2) = X.z();
        P(2 * i, 3) = 1;
        P(2 * i, 4) = 0;
        P(2 * i, 5) = 0;
        P(2 * i, 6) = 0;
        P(2 * i, 7) = 0;
        P(2 * i, 8) = -x.x() * X.x();
        P(2 * i, 9) = -x.x() * X.y();
        P(2 * i, 10) = -x.x() * X.z();
        P(2 * i, 11) = -x.x();

        // Y-coordinate row
        P(2 * i + 1, 0) = 0;
        P(2 * i + 1, 1) = 0;
        P(2 * i + 1, 2) = 0;
        P(2 * i + 1, 3) = 0;
        P(2 * i + 1, 4) = X.x();
        P(2 * i + 1, 5) = X.y();
        P(2 * i + 1, 6) = X.z();
        P(2 * i + 1, 7) = 1;
        P(2 * i + 1, 8) = -x.y() * X.x();
        P(2 * i + 1, 9) = -x.y() * X.y();
        P(2 * i + 1, 10) = -x.y() * X.z();
        P(2 * i + 1, 11) = -x.y();
    }

    // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.
    Matrix U, D, V;
    svd_decompose(P, U, D, V);

    // Extract m from the last column of V
    Vector m = V.get_column(V.cols() - 1);

    // Reformat vector m into matrix M
    Matrix33 M;
    M.set_column(0, {m[0], m[1], m[2]});
    M.set_column(1, {m[4], m[5], m[6]});
    M.set_column(2, {m[8], m[9], m[10]});
    t = Vector3D(m[3], m[7], m[11]); // Translation vector

    // TODO: extract intrinsic parameters from M.
    Vector3D a1 = M.get_row(0);
    Vector3D a2 = M.get_row(1);
    Vector3D a3 = M.get_row(2);
    double rho = 1 / norm(a3);
    cx = rho * rho * dot(a1, a3);
    cy = rho * rho * dot(a2, a3);
    double cosTheta = -dot(cross(a1, a3), cross(a2, a3)) / (norm(cross(a1, a3)) * norm(cross(a2, a3)));
    double sinTheta = sqrt(1 - cosTheta * cosTheta);
    fx = rho * rho * norm(cross(a1, a3)) * sinTheta;
    fy = rho * rho * norm(cross(a2, a3)) * sinTheta;

    // Calculate extrinsic parameters
    Vector3D r1 = cross(a2, a3) / norm(cross(a2, a3));
    Vector3D r3 = rho * a3;
    Vector3D r2 = cross(r3, r1);

    R.set_column(0, r1);
    R.set_column(1, r2);
    R.set_column(2, r3);

    std::cout << "Calibration successful, updating camera parameters.\n";

    std::cout << "\n\tTODO: After you implement this function, please return 'true' - this will trigger the viewer to\n"
                 "\t\tupdate the rendering using your recovered camera parameters. This can help you to visually check\n"
                 "\t\tif your calibration is successful or not.\n\n" << std::flush;
    return true;
}

















