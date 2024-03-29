# -*- coding: utf-8 -*
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

from BoltzTraP2.sphere.frontend import calc_reciprocal_stars
from BoltzTraP2.sphere.frontend import calc_sphere_quotient_set
from BoltzTraP2.sphere.frontend import calc_reciprocal_degeneracies
from BoltzTraP2.misc import info

import math
import numpy as np
import numpy.linalg as la

import spglib


def compute_bounds(lattvec, r):
    """Obtain bounds on each of the three components of a column vector
    J = [j1, j2, j3]  such that |lattvec @ J| <= r.
    """
    if lattvec.shape != (3, 3):
        raise ValueError("lattvec must be a 3 x 3 matrix")
    try:
        la.inv(lattvec)
    except la.LinAlgError:
        raise ValueError("singular lattice matrix")
    if r <= 0:
        raise ValueError("the radius must be positive")
    nruter = np.zeros(3, dtype=int)
    # We look for a set of lattice planes with only one non-zero index
    # that can completely contain the sphere. It is evident that those
    # always exist, and in fact we only need to find an arbitrary set of planes
    # that fall outside the sphere.
    # We start by looking for the planes tangent to the sphere, i.e, those whose
    # minimum distance to the origin is equal to the radius.
    # The minimum distance of a plane to the origin is the modulus of the
    # projection of any point of the plane over its unit normal, i.e.,
    # |n @ r| / |n|. In this particular case, r is a multiple of a single
    # lattice vector, and the normal n is the vector from reciprocal space
    # with the same reduced coordinates. With that in mind, what follows
    # is trivially derived from the condition
    # |n @ r| ** 2  / |n| ** 2 == r ** 2.
    invmetric = la.inv(lattvec.T @ lattvec)
    for k in range(3):
        u = np.eye(3)[:, k]
        tangent = r * np.sqrt(u @ invmetric @ u)
        # We round up the solution to find an integer bound.
        nruter[k] = math.ceil(tangent)
    return nruter


def compute_radius(atoms, nkpt, symprec=1e-4):
    """Estimate the right radius to be passed to Equivalence_builder to
    obtain the desired number of points.
    """
    lattvec = atoms.get_cell().T
    vol = abs(np.linalg.det(lattvec))
    spg = spglib.get_symmetry_dataset(atoms, symprec)
    rotations = spg["rotations"].astype(int)
    nrot = len(rotations)
    info("{} rotations".format(nrot))
    newrotations = set()
    # Include time-reversal symmetry in the estimate if total inversion of
    # the axes is not among the symmetries
    inv = 2
    for i in range(nrot):
        tmp = tuple(tuple(j) for j in rotations[i].tolist())
        newrotations.add(tmp)
        if tmp == ((-1, 0, 0), (0, -1, 0), (0, 0, -1)):
            inv = 1
    nrot = len(newrotations)
    info(nrot, "unique rotations", "including the inversion"
         if inv == 1 else "")
    npoints = nkpt * nrot * inv
    info(npoints, "total k points in the sphere")
    radius = (3. / (4. * np.pi) * npoints * vol)**(1. / 3.)
    return radius


def get_equivalences(atoms, nkpt):
    """Get a list of approximately nkpt equivalence classes of k points for the
    provided atoms object.
    """
    radius = compute_radius(atoms, nkpt)
    bounds = compute_bounds(atoms.get_cell().T, radius)
    return calc_sphere_quotient_set(atoms, radius, bounds)
