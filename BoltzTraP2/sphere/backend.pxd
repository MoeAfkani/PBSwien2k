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

# Wrappers for the C++ code in backend.cpp
from libcpp.vector cimport vector

cdef extern from "backend.hpp":
    cdef cppclass Sphere_equivalence_builder:
        Sphere_equivalence_builder(double[3][3],
                                   double[][3],
                                   int[],
                                   int,
                                   double,
                                   int[3],
                                   double) except +
        vector[int] get_mapping()
        void get_point(int, int[3])
    cdef cppclass Degeneracy_counter:
        Degeneracy_counter(double[3][3],
                           double[][3],
                           int[],
                           int,
                           double) except +
        vector[int] get_mapping()
        int count_degeneracy(double[3])
        void get_point(int, double[3])
