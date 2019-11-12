//    BoltzTraP2, a program for interpolating band structures and
//    calculating
//                semi-classical transport coefficients.
//    Copyright (C) 2017 Georg K. H. Madsen <georg.madsen@tuwien.ac.at>
//    Copyright (C) 2017 Jesús Carrete
//    <jesus.carrete.montana@tuwien.ac.at>
//    Copyright (C) 2017 Matthieu J. Verstraete
//    <matthieu.verstraete@ulg.ac.be>
//
//    This file is part of BoltzTraP2.
//
//    BoltzTraP2 is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published
//    by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    BoltzTraP2 is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with BoltzTraP2.  If not, see
//    <http://www.gnu.org/licenses/>.

#include "backend.hpp"

#include <cmath>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <set>
#include <utility>

extern "C" {
#include <spglib.h>
}

/// Premultiply a 3x3 matrix by its transpose
///
/// @param[in] m - a 3x3 matrix
/// @param[out] dest - matrix to store the result
template<typename T> void calc_square_3m(const T m[3][3], T dest[3][3])
{
    for (int i = 0; i < 3; ++i) {
        for (int j = i; j < 3; ++j) {
            dest[i][j] = 0;
        }
    }
    for (int l = 0; l < 3; ++l) {
        for (int i = 0; i < 3; ++i) {
            for (int j = i; j < 3; ++j) {
                dest[i][j] += m[l][i] * m[l][j];
            }
        }
    }
    for (int i = 1; i < 3; ++i) {
        for (int j = 0; j < i; ++j) {
            dest[i][j] = dest[j][i];
        }
    }
}

/// Multiply a 3x3 matrix by a 3-vector
///
/// @param[in] m - a 3x3 matrix
/// @param[in] v - a 3-vector
/// @param[out] dest - array to store the result
template<typename T1, typename T2>
void calc_3mv(const matrix33<T1>& m, const vector3<T2>& v,
              vector3<T2>& dest)
{
    for (int i = 0; i < 3; ++i) {
        dest[i] = 0;
        for (int j = 0; j < 3; ++j) {
            dest[i] += m[i][j] * v[j];
        }
    }
}

/// Multiply two 3x3 matrices
///
/// @param[in] x - first matrix
/// @param[in] y - second matrix
/// @param[out] dest - array to store the result
template<typename T1, typename T2>
void calc_3mm(const T1 x[3][3], const matrix33<T2>& y,
              matrix33<T1>& dest)
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            dest[i][j] = 0.;
            for (int k = 0; k < 3; ++k) {
                dest[i][j] += x[i][k] * y[k][j];
            }
        }
    }
}

/// Multiply two 3x3 matrices
///
/// @param[in] x - first matrix
/// @param[in] y - second matrix
/// @param[out] dest - array to store the result
template<typename T1, typename T2>
void calc_3mm(const matrix33<T1>& x, const T2 y[3][3],
              matrix33<T1>& dest)
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            dest[i][j] = 0.;
            for (int k = 0; k < 3; ++k) {
                dest[i][j] += x[i][k] * y[k][j];
            }
        }
    }
}

/// Multiply two 3x3 matrices
///
/// @param[in] x - first matrix
/// @param[in] y - second matrix
/// @param[out] dest - array to store the result
template<typename T1, typename T2>
void calc_3mm(const matrix33<T1>& x, const matrix33<T2>& y,
              matrix33<T1>& dest)
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            dest[i][j] = 0.;
            for (int k = 0; k < 3; ++k) {
                dest[i][j] += x[i][k] * y[k][j];
            }
        }
    }
}

/// Obtain the transposed inverse of a 3x3 matrix
///
/// @param[in] m - a 3x3 matrix
/// @param[out] dest - matriix to store the result
template<typename T>
void calc_invt_3m(const T m[3][3], matrix33<T>& dest)
{
    T det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
            m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
            m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    dest[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) / det;
    dest[1][0] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) / det;
    dest[2][0] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det;
    dest[0][1] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) / det;
    dest[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det;
    dest[2][1] = -(m[0][0] * m[1][2] - m[1][0] * m[0][2]) / det;
    dest[0][2] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) / det;
    dest[1][2] = -(m[0][0] * m[2][1] - m[2][0] * m[0][1]) / det;
    dest[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) / det;
}


/// Transpose a 3x3 matrix
///
/// @param[in] m - a 3x3 matrix
/// @param[out] dest - matrix to store the result
template<typename T> void calc_t_3m(const matrix33<T>& m, T dest[3][3])
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            dest[i][j] = m[j][i];
        }
    }
}

/// Transpose a 3x3 matrix
///
/// @param[in] m - a 3x3 matrix
/// @param[out] dest - matrix to store the result
template<typename T>
void calc_t_3m(const matrix33<T>& m, matrix33<T>& dest)
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            dest[i][j] = m[j][i];
        }
    }
}

/// Compare two 3-vectors for strict equality
///
/// @param[in] v1 - first vector
/// @param[in] v2 - second vector
/// @return true if the two vectors are equal or false otherwise
template<typename T>
bool match_3m(const vector3<T>& v1, const vector3<T>& v2)
{
    return v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2];
}

/// Compare two 3-vectors for modular equality within a tolerance
///
/// @param[in] v1 - first vector
/// @param[in] v2 - second vector
/// @param[in] epsilon - tolerance
/// @return true if the difference between both vectors is a vector of
/// integers up to the tolerance, or false otherwise
template<typename T>
bool match_3m_reciprocal(const vector3<T>& v1, const vector3<T>& v2,
                         double epsilon)
{
    vector3<T> delta({v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]});
    return std::abs(delta[0] - std::round(delta[0])) < epsilon &&
           std::abs(delta[1] - std::round(delta[1])) < epsilon &&
           std::abs(delta[2] - std::round(delta[2])) < epsilon;
}

Symmetry_analyzer::Symmetry_analyzer(double lattvec[3][3],
                                     double positions[][3],
                                     const int types[], int natoms,
                                     double symprec)
{
    // Obtain the rotation matrices
    this->analyze_symmetries(lattvec, positions, types, natoms,
                             symprec);
}

void Symmetry_analyzer::analyze_symmetries(double lattvec[3][3],
                                           double positions[][3],
                                           const int types[],
                                           int natoms, double symprec)
{
    // Call spglib to get all information about symmetries
    SpglibDataset* data =
        spg_get_dataset(lattvec, positions, types, natoms, symprec);
    if (data == nullptr) {
        throw Sphere_exception("spglib's spg_get_dataset"
                               " returned a NULL pointer");
    }
    // Scan through the operations and keep only a unique set of
    // rotations plus the corresponding roto-inversions.
    auto nops = data->n_operations;
    std::set<matrix33<int>> unique;
    for (int i = 0; i < nops; ++i) {
        matrix33<int> tmp;
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                tmp[r][c] = data->rotations[i][r][c];
            }
        }
        unique.emplace(tmp);
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                tmp[r][c] = -tmp[r][c];
            }
        }
        unique.emplace(tmp);
    }
    this->rotations.assign(unique.begin(), unique.end());
    // Free up the space allocated by spglib
    spg_free_dataset(data);
}

/// Simple class for computing squared norms of lattice vectors
class Lattice_tape
{
public:
    /// Basic constructor.
    ///
    /// @param[in] lattvec - lattice vectors, as columns
    Lattice_tape(const double lattvec[3][3])
    {
        calc_square_3m(lattvec, this->metric);
    }
    /// Compute the squared norm of a lattice vector
    ///
    /// @param[in] v - vector in direct coordinates
    /// @return the squared norm of v
    double measure(const vector3<int>& v) const
    {
        double nruter = 0.;
        double dv[3];
        for (int i = 0; i < 3; ++i) {
            dv[i] = static_cast<double>(v[i]);
        }
        for (int i = 0; i < 3; ++i) {
            nruter += this->metric[i][i] * dv[i] * dv[i];
            for (int j = i + 1; j < 3; ++j) {
                nruter += 2. * this->metric[i][j] * dv[i] * dv[j];
            }
        }
        return nruter;
    }

private:
    /// Metric tensor
    double metric[3][3];
};

/// Comparison function for pairs or tuples based on the second element
///
/// @param p1 - first operand
/// @param p2 - second operand
/// @return the result of std::get<1>(p1) < std::get<1>(p2)
template<typename T> bool compare_second(const T& p1, const T& p2)
{
    return std::get<1>(p1) < std::get<1>(p2);
}

Sphere_equivalence_builder::Sphere_equivalence_builder(
    double lattvec[3][3], double positions[][3], const int types[],
    int natoms, double r, const int bounds[3], double symprec)
    : Equivalence_builder(lattvec, positions, types, natoms, symprec)
{
    // Generate a list of all lattice points in the intersection and
    // store their squared norms for later use
    this->create_grid(lattvec, r, bounds);
    // Initialize the point mapping
    this->mapping.resize(this->grid.size());
    std::fill(this->mapping.begin(), this->mapping.end(), -1);
    // Find an equivalence class for each point
    auto lo = std::make_pair(vector3<int>({0, 0, 0}), 0.);
    auto hi = std::make_pair(vector3<int>({0, 0, 0}), 0.);
    vector3<int> image;
    for (std::size_t i = 0; i < this->grid.size(); ++i) {
        // Each point can only be equivalent to a point with
        // the same norm. Some room is left for rounding errors
        auto& point = std::get<0>(this->grid[i]);
        auto sqnorm = std::get<1>(this->grid[i]);
        std::get<1>(lo) = (1. - symprec) * sqnorm;
        std::get<1>(hi) = (1. + symprec) * sqnorm;
        std::size_t lbound = std::distance(
            this->grid.begin(),
            std::lower_bound(this->grid.begin(), this->grid.end(), lo,
                             compare_second<point_with_sqnorm>));
        std::size_t ubound = std::distance(
            this->grid.begin(),
            std::upper_bound(this->grid.begin(), this->grid.end(), hi,
                             compare_second<point_with_sqnorm>));
        ubound = std::min(ubound, i);
        for (std::size_t o = 0; o < this->rotations.size(); ++o) {
            // Apply each symmetry operation to the point
            calc_3mv(this->rotations[o], point, image);
            // And compare the result to each of the candidates
            for (std::size_t j = lbound; j < ubound; ++j) {
                if (match_3m(image, std::get<0>(this->grid[j]))) {
                    // If they match, find the representative of the
                    // equivalence class, assign this point to the class
                    // and exit the loop
                    int k = j;
                    while (this->mapping[k] != k) {
                        k = this->mapping[k];
                    }
                    this->mapping[i] = k;
                    break;
                }
            }
            if (this->mapping[i] != -1) {
                break;
            }
        }
        // If no equivalence class was found for this point, start a new
        // one
        if (this->mapping[i] == -1) {
            this->mapping[i] = i;
        }
    }
}

void Sphere_equivalence_builder::create_grid(double lattvec[3][3],
                                             double r,
                                             const int bounds[3])
{
    double r2(r * r);
    Lattice_tape tape(lattvec);
    for (int i = -bounds[0]; i <= bounds[0]; ++i) {
        for (int j = -bounds[1]; j <= bounds[1]; ++j) {
            for (int k = -bounds[2]; k <= bounds[2]; ++k) {
                vector3<int> point({i, j, k});
                double n2 = tape.measure(point);
                if (n2 < r2) {
                    this->grid.emplace_back(std::make_pair(point, n2));
                }
            }
        }
    }
    // The list is stored in order of increasing norm
    std::sort(this->grid.begin(), this->grid.end(),
              compare_second<point_with_sqnorm>);
}

Degeneracy_counter::Degeneracy_counter(double lattvec[3][3],
                                       double positions[][3],
                                       const int types[], int natoms,
                                       double symprec)
    : Equivalence_builder(lattvec, positions, types, natoms, symprec),
      tolerance(symprec)
{
    // Obtain the rotation matrices in the reciprocal basis
    double metric[3][3];
    calc_square_3m(lattvec, metric);
    matrix33<double> itmetric;
    calc_invt_3m(metric, itmetric);
    double imetric[3][3];
    calc_t_3m(itmetric, imetric);
    for (const auto& r : this->rotations) {
        matrix33<double> tmp;
        calc_3mm(metric, r, tmp);
        matrix33<double> kr;
        calc_3mm(tmp, imetric, kr);
        this->krotations.emplace_back(kr);
    }
}

int Degeneracy_counter::count_degeneracy(double point[3])
{
    // Fill the list with all possible images of the point.
    this->kpoints.clear();
    vector3<double> kpoint({point[0], point[1], point[2]});
    vector3<double> image;
    for (std::size_t o = 0; o < this->krotations.size(); ++o) {
        calc_3mv(this->krotations[o], kpoint, image);
        for (int j = 0; j < 3; ++j) {
            image[j] -= std::round(image[j]);
        }
        this->kpoints.emplace_back(image);
    }
    // Initialize the point mapping
    this->mapping.resize(this->kpoints.size());
    std::fill(this->mapping.begin(), this->mapping.end(), -1);
    // Build the mapping by finding points related by translations
    for (std::size_t i = 0; i < this->kpoints.size(); ++i) {
        auto& point = this->kpoints[i];
        for (std::size_t j = 0; j < i; ++j) {
            if (match_3m_reciprocal(point, this->kpoints[j],
                                    this->tolerance)) {
                int k = j;
                while (this->mapping[k] != k) {
                    k = this->mapping[k];
                }
                this->mapping[i] = k;
                break;
            }
        }
        if (this->mapping[i] == -1) {
            this->mapping[i] = i;
        }
    }
    /// Count and return the number of classes in the mapping
    std::sort(this->mapping.begin(), this->mapping.end());
    return std::unique(this->mapping.begin(), this->mapping.end()) -
           this->mapping.begin();
}

