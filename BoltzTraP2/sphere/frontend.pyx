#    BoltzTraP2, a program for interpolating band structures and calculating
#                semi-classical transport coefficients.
#    Copyright (C) 2017 Georg K. H. Madsen <georg.madsen@tuwien.ac.at>
#    Copyright (C) 2017 Jesús Carrete <jesus.carrete.montana@tuwien.ac.at>
#    Copyright (C) 2017 Matthieu J. Verstraete <matthieu.verstraete@ulg.ac.be>
#
#    This file is part of BoltzTraP2.
#
#    BoltzTraP2 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    BoltzTraP2 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with BoltzTraP2.  If not, see <http://www.gnu.org/licenses/>.

# cython: boundscheck=False
# cython: cdivision=True

from BoltzTraP2.sphere.backend cimport Sphere_equivalence_builder
from BoltzTraP2.sphere.backend cimport Degeneracy_counter

import numpy as np
import numpy.linalg as la

from libc.stdlib cimport malloc
from libc.stdlib cimport free
from libcpp.vector cimport vector
cimport numpy as np

np.import_array()


def calc_sphere_quotient_set(atoms, radius, bounds, symprec=1e-4):
    """
    Return a set of equivalence classes for the lattice points in the
    intersection between the supercell based on the lattice of "atoms" with
    the maximum absolute coordinates determined by "bounds", and a sphere of
    radius "radius", with respect to the rotations of the structure in "atoms"
    plus an inversion about the origin.
    """
    lattvec = atoms.get_cell().T
    # Translate the information in the atoms object into a set of
    # structures compatible with spglib's C API
    cdef int natoms = len(atoms.numbers)
    cdef int i
    cdef int j
    cdef double c_lattvec[3][3]
    for i in range(3):
        for j in range(3):
            c_lattvec[i][j] = lattvec[i, j]
    positions = atoms.get_scaled_positions()
    cdef double(*c_positions)[3]
    c_positions = <double(*)[3] > malloc(natoms * sizeof(double[3]))
    if c_positions is NULL:
        raise MemoryError
    for i in range(natoms):
        for j in range(3):
            c_positions[i][j] = positions[i][j]
    cdef int[:] c_types = np.ascontiguousarray(atoms.numbers, dtype=np.intc)
    cdef int[:] c_bounds = np.ascontiguousarray(bounds, dtype=np.intc)
    # Pass these to Sphere_equivalence_builder, which does the actual work
    eq_builder = new Sphere_equivalence_builder(c_lattvec,
                                                c_positions,
                                                & (c_types[0]),
                                                natoms,
                                                radius,
                                                & (c_bounds[0]),
                                                symprec)
    free(c_positions)
    # Get a mapping from the whole point set to a set of irreducible points
    mapping = np.array(eq_builder.get_mapping())
    # And the coordinates of the points in the grid
    grid = np.empty((mapping.size, 3), order="C", dtype=np.intc)
    cdef int[:, :] c_grid = grid
    for i in range(mapping.size):
        eq_builder.get_point(i, & (c_grid[i, 0]))
    # Build the quotient group with the survivors as representatives.
    # Step 1: put equivalent points together.
    indices = mapping.argsort()
    grid = grid[indices, :]
    mapping = mapping[indices]
    # Step 2: find the indices where the mapping changes.
    splitpoints = np.nonzero((mapping[:-1] != mapping[1:]))[0] + 1
    # Step 3: apply the scissors there.
    equivalences = np.split(grid, splitpoints, axis=0)
    return equivalences


def calc_reciprocal_degeneracies(atoms, kpoints, symprec=1e-4):
    """
    For each reciprocal-space point in the nx3 array kpoints, count the
    number of equivalent reciprocal-space points in a single copy of the
    Brillouin zone. Return an array with those counts.
    """
    lattvec = atoms.get_cell().T
    # Translate the information in the atoms object into a set of
    # structures compatible with spglib's C API
    cdef int natoms = len(atoms.numbers)
    cdef int npoints = kpoints.shape[0]
    cdef int i
    cdef int j
    cdef double c_lattvec[3][3]
    for i in range(3):
        for j in range(3):
            c_lattvec[i][j] = lattvec[i, j]
    positions = atoms.get_scaled_positions()
    cdef double(*c_positions)[3]
    c_positions = <double(*)[3] > malloc(natoms * sizeof(double[3]))
    if c_positions is NULL:
        raise MemoryError
    for i in range(natoms):
        for j in range(3):
            c_positions[i][j] = positions[i][j]
    cdef int[:] c_types = np.ascontiguousarray(atoms.numbers, dtype=np.intc)
    # Pass these to Degeneracy_counter, which does the actual work
    eq_builder = new Degeneracy_counter(c_lattvec,
                                        c_positions,
                                        & (c_types[0]),
                                        natoms,
                                        symprec)
    free(c_positions)
    # Interrogate the object about the degeneracy of each k point and
    # return the result.
    cdef double[:, :] c_kpoints = np.ascontiguousarray(kpoints)
    nruter = np.empty(npoints, dtype=np.intc)
    for i in range(npoints):
        nruter[i] = eq_builder.count_degeneracy(& (c_kpoints[i, 0]))
    return nruter


def calc_reciprocal_stars(atoms, kpoints, symprec=1e-4):
    """
    For each reciprocal-space point in the nx3 array kpoints, find all the
    equivalent reciprocal-space points in the first copy of the Brillouin zone.
    Return a list of arrays with those points.
    """
    lattvec = atoms.get_cell().T
    # Translate the information in the atoms object into a set of
    # structures compatible with spglib's C API
    cdef int natoms = len(atoms.numbers)
    cdef int npoints = kpoints.shape[0]
    cdef int i
    cdef int j
    cdef double c_lattvec[3][3]
    for i in range(3):
        for j in range(3):
            c_lattvec[i][j] = lattvec[i, j]
    positions = atoms.get_scaled_positions()
    cdef double(*c_positions)[3]
    c_positions = <double(*)[3] > malloc(natoms * sizeof(double[3]))
    if c_positions is NULL:
        raise MemoryError
    for i in range(natoms):
        for j in range(3):
            c_positions[i][j] = positions[i][j]
    cdef int[:] c_types = np.ascontiguousarray(atoms.numbers, dtype=np.intc)
    # Pass these to Degeneracy_counter, which does the actual work
    eq_builder = new Degeneracy_counter(c_lattvec,
                                        c_positions,
                                        & (c_types[0]),
                                        natoms,
                                        symprec)
    free(c_positions)
    # Interrogate the object about the degeneracy of each k point and
    # find a representative of each class in the generated mapping.
    # Each of these classes is a set of points related by translations
    # in reciprocal space by a reciprocal-lattice vector.
    cdef double[:, :] c_kpoints = np.ascontiguousarray(kpoints)
    nruter = []
    cdef double[:, :] c_reprs
    for i in range(npoints):
        eq_builder.count_degeneracy(& (c_kpoints[i, 0]))
        mapping = np.array(eq_builder.get_mapping())
        positions = np.unique(mapping)
        reprs = np.empty((positions.size, 3), order="C")
        c_reprs = reprs
        for j, p in enumerate(positions):
            eq_builder.get_point(p, & (c_reprs[j, 0]))
        nruter.append(reprs)
    return nruter
