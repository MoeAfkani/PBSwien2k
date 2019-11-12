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

#pragma once

#include <array>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

/// Basic exception class used for all errors in this module
class Sphere_exception : public std::runtime_error
{
public:
    /// Basic constructor
    Sphere_exception(const std::string& message)
        : std::runtime_error(message)
    {
    }
};

/// Type declaration for basic 3-vectors for internal use
template<typename T> using vector3 = std::array<T, 3>;

/// Type declaration for basic 3x3 matrices for internal use
template<typename T> using matrix33 = std::array<vector3<T>, 3>;

/// Type representing a lattice point with its squared norm attached
using point_with_sqnorm = std::pair<vector3<int>, double>;

// Base class for classes that deal with symmetries in this module.
class Symmetry_analyzer
{
public:
    /// Basic constructor
    ///
    /// @param[in] lattvec - lattice vectors, as columns
    /// @param[in] positions - atomic positions in direct coordinates
    /// @param[in] types - list of atom types
    /// @param[in] natoms - number of atoms in the structure
    /// @param[in] symprec - tolerance for symmetry search
    Symmetry_analyzer(double lattvec[3][3], double positions[][3],
                      const int types[], int natoms, double symprec);

protected:
    /// Set of rotations from the space group
    std::vector<matrix33<int>> rotations;
    /// Analyze the structure with spglib and store the rotations
    /// of the space group
    ///
    /// @param[in] lattvec - lattice vectors, as columns
    /// @param[in] positions - atomic positions in direct coordinates
    /// @param[in] types - list of atom types
    /// @param[in] natoms - number of atoms in the structure
    /// @param[in] symprec - tolerance for symmetry search
    void analyze_symmetries(double lattvec[3][3], double positions[][3],
                            const int types[], int natoms,
                            double symprec);
};

/// Base class for classes that build equivalences between points
class Equivalence_builder : public Symmetry_analyzer
{
public:
    /// Basic constructor
    ///
    /// @param[in] lattvec - lattice vectors, as columns
    /// @param[in] positions - atomic positions in direct coordinates
    /// @param[in] types - list of atom types
    /// @param[in] natoms - number of atoms in the structure
    /// @param[in] symprec - tolerance for symmetry search
    Equivalence_builder(double lattvec[3][3], double positions[][3],
                        const int types[], int natoms, double symprec)
        : Symmetry_analyzer(lattvec, positions, types, natoms, symprec)
    {
    }
    /// Return a copy of the point mapping
    std::vector<int> get_mapping() const
    {
        return this->mapping;
    }

protected:
    /// Index of the equivalence class of each point of the grid
    std::vector<int> mapping;
};

/// Class used to classify equivalent lattice points in the intersection
/// between a supercell and a sphere
class Sphere_equivalence_builder : public Equivalence_builder
{
public:
    /// Basic constructor
    ///
    /// @param[in] lattvec - lattice vectors, as columns
    /// @param[in] positions - atomic positions in direct coordinates
    /// @param[in] types - list of atom types
    /// @param[in] natoms - number of atoms in the structure
    /// @param[in] r - radius of the sphere
    /// @param[in] bounds - maximum absolute value of each direct
    /// coordinate
    /// in the supercell
    /// @param[in] symprec - tolerance for symmetry search
    Sphere_equivalence_builder(double lattvec[3][3],
                               double positions[][3], const int types[],
                               int natoms, double r,
                               const int bounds[3], double symprec);
    /// Load the coordinates of a grid point into the argument
    ///
    /// Beware that the index is not checked
    /// @param[in] index - index of the point in the grid
    /// @param[out] point - array to store the coordinates of the point
    void get_point(int index, int point[3]) const
    {
        auto& gpoint = std::get<0>(this->grid[index]);
        for (int i = 0; i < 3; ++i) {
            point[i] = gpoint[i];
        }
    }

private:
    /// List of lattice points with their squared norm attached
    std::vector<point_with_sqnorm> grid;
    /// Fill the "grid" vector with all points in the intersection of
    /// the
    /// supercell and the sphere
    ///
    /// @param[in] lattvec - lattice vectors, as columns
    /// @param[in] r - radius of the sphere
    /// @param[in] bounds - maximum absolute value of each direct
    /// coordinate in the supercell
    void create_grid(double lattvec[3][3], double r,
                     const int bounds[3]);
};

/// Class used to find all points in the BZ related
/// to a single point by symmetry.
class Degeneracy_counter : public Equivalence_builder
{
public:
    /// Basic constructor
    ///
    /// @param[in] lattvec - lattice vectors, as columns
    /// @param[in] positions - atomic positions in direct coordinates
    /// @param[in] types - list of atom types
    /// @param[in] natoms - number of atoms in the structure
    /// @param[in] symprec - tolerance for symmetry search
    Degeneracy_counter(double lattvec[3][3], double positions[][3],
                       const int types[], int natoms, double symprec);
    /// Obtain the number of symmetric versions of a point in reciprocal
    /// space, neglecting simple translations by a reciprocal lattice
    /// vector
    ///
    /// @param[in] point - point to analyze
    /// @return the degeneracy of the point
    int count_degeneracy(double point[3]);
    /// Load the coordinates of a point into the argument
    ///
    /// Beware that the index is not checked
    /// @param[in] index - index of the point in the grid
    /// @param[out] point - array to store the coordinates of the point
    void get_point(int index, double point[3]) const
    {
        for (int i = 0; i < 3; ++i) {
            point[i] = this->kpoints[index][i];
        }
    }

private:
    /// Tolerance for symmetry search
    double tolerance;
    /// Set of rotations from the space group, for vectors in reciprocal
    /// coordinates
    std::vector<matrix33<double>> krotations;
    /// List of reciprocal-space points
    std::vector<vector3<double>> kpoints;
};
