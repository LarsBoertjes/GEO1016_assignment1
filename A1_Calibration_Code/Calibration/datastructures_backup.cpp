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

    /// Below are a few examples showing some useful data structures and functions.

    // This is a 1D array of 'double' values. Alternatively, you can use 'double mat[25]' but you cannot change it
    // length. With 'std::vector', you can append/delete/insert elements, and much more. The 'std::vector' can store
    // not only 'double', but also any other types of objects. In case you may want to learn more about 'std::vector'
    // check here: https://en.cppreference.com/w/cpp/container/vector
    //std::vector<double> array = {1, 3, 3, 4, 7, 6, 2, 8, 2, 8, 3, 2, 4, 9, 1, 7, 3, 23, 2, 3, 5, 2, 1, 5, 8, 9, 22};
    //array.push_back(5); // append 5 to the array (so the size will increase by 1).
    //array.insert(array.end(), 10, 3);  // append ten 3 (so the size will grow by 10).

    /// To access the value of an element.
    //double a = array[2];

    /// define a 2D vector/point
    //Vector2D b(1.1, 2.2);

    /// define a 3D vector/point
    //Vector3D c(1.1, 2.2, 3.3);

    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    //Vector2D p = c.cartesian();

    /// get the Homogeneous coordinates of p
    //Vector3D q = p.homogeneous();

    /// the length of a vector
    //double len = p.length();
    /// the squared length of a vector
    //double sqr_len = p.length2();

    /// the dot product of two vectors
    //double dot_prod = dot(p, q);

    /// the cross product of two vectors
    //Vector cross_prod = cross(c, q);

    /// normalize this vector
    //cross_prod.normalize();

    // Define an m-by-n double valued matrix.
    // Here I use the above array to initialize it. You can also use A(i, j) to initialize/modify/access its elements.
    //const int m = 6, n = 5;
    //Matrix A(m, n, array.data());    // 'array.data()' returns a pointer to the array.
//    std::cout << "M: \n" << A << std::endl;

    /// define a 3 by 4 matrix (and all elements initialized to 0.0)
    //Matrix M(3, 4, 0.0);

    /// set first row by a vector
    //M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
    //M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    //Matrix33 B;

    /// define and initialize a 3 by 3 matrix
    /*Matrix33 T(1.1, 2.2, 3.3,
               0, 2.2, 3.3,
               0, 0, 1);*/

    /// define and initialize a 3 by 4 matrix
    /*Matrix34 P(1.1, 2.2, 3.3, 0,
               0, 2.2, 3.3, 1,
               0, 0, 1, 1);*/

    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    //Matrix W(15, 9, 0.0);
    /// set the first row by a 9-dimensional vector
    //W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>

    /// get the number of rows.
    //int num_rows = W.rows();

    /// get the number of columns.
    //int num_cols = W.cols();

    /// get the the element at row 1 and column 2
    //double value = W(1, 2);

    /// get the last column of a matrix
    //Vector last_column = W.get_column(W.cols() - 1);

    /// define a 3 by 3 identity matrix
    //Matrix33 I = Matrix::identity(3, 3, 1.0);

    /// matrix-vector product
    //Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

    //Matrix U(m, m, 0.0);   // initialized with 0s
    //Matrix S(m, n, 0.0);   // initialized with 0s
    //Matrix V(n, n, 0.0);   // initialized with 0s

    // Compute the SVD decomposition of A
    //svd_decompose(A, U, S, V);

    // Now let's check if the SVD result is correct

    // Check 1: U is orthogonal, so U * U^T must be identity
//    std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;

    // Check 2: V is orthogonal, so V * V^T must be identity
//    std::cout << "V*V^T: \n" << V * transpose(V) << std::endl;

    // Check 3: S must be a diagonal matrix
//    std::cout << "S: \n" << S << std::endl;

    // Check 4: according to the definition, A = U * S * V^T
//    std::cout << "M - U * S * V^T: \n" << A - U * S * transpose(V) << std::endl;

    // Compute the inverse of a matrix
    //Matrix invT;
    //inverse(T, invT);
    // Let's check if the inverse is correct
//    std::cout << "B * invB: \n" << B * invB << std::endl;

    // TODO: the above code just demonstrates some useful data structures and APIs. Please remove all above code in your
    //       final submission.

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    std::cout << "\n[Liangliang]:\n"
                 "\tThe input parameters of this function are:\n"
                 "\t\t- points_3d: An array of 3D points (input to this function)\n"
                 "\t\t- points_2d: An array of 2D image points (input to this function)\n"
                 "\tThis function must return either 'true' on success or 'false' otherwise. On success, the camera\n"
                 "\tparameters are returned by the following variables:\n"
                 "\t\t- fx and fy: the focal lengths\n"
                 "\t\t- cx and cy: the principal point\n"
                 "\t\t- s: the skew factor, i.e., s = -alpha * cot(theta)\n"
                 "\t\t- R: the 3x3 rotation matrix encoding camera orientation\n"
                 "\t\t- t: a 3D vector encoding camera location.\n"
                 "\tIMPORTANT: don't forget to write your recovered parameters to the above variables." << std::endl;

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
    Matrix U, S, V;
    svd_decompose(P, U, S, V);

    // Extract m from the last column of V
    Vector m = V.get_column(V.cols() - 1);
    std::cout << m.size() << std::endl;

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

















