# -*- coding: utf-8 -*-
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

import xml.etree.ElementTree as et

import numpy as np
import scipy as sp
import scipy.spatial.distance
import ase

from BoltzTraP2.misc import ffloat
from BoltzTraP2.units import *


def read_GENE_eneandmat(filename):
    """Read in a .energy file in GENE format, describing the electronic bands.

    The file contains information about k points, band energies, and
    (optionally) momentum matrix elements.

    Args:
        filename: path to the .energy file

    Returns:
        A 5-tuple. The first element is the Fermi level. The second is the
        maximum occupancy per state. The third is an (nkpoints, 3) array with
        the coordinates of the irreducible k points of the system. The fourth
        contains the energy bands in an array. The fifth is either an array
        with the derivatives of the bands, or None if such information is not
        contained in the file.
    """
    lines = open(filename, "r").readlines()
    # Line 1: title string
    # Line 2: nk, nspin, Fermi level(Ry)
    linenumber = 1
    tmp = lines[linenumber].split()
    nk, nspin, efermi = int(tmp[0]), int(tmp[1]), ffloat(tmp[2])

    minband = np.infty
    ebands1 = []
    mommat1 = []
    kpoints = []
    for ispin in range(nspin):
        for ik in range(nk):
            # k block: line 1 = kx ky kz nband
            linenumber += 1
            tmp = lines[linenumber].split()
            nband = int(tmp[3])
            if nband < minband:
                minband = nband
            kpoints.append([ffloat(i) for i in tmp[0:3]])
            eband = []
            vband = []
            for ib in range(nband):
                linenumber += 1
                fields = lines[linenumber].split()
                e = ffloat(fields[0])
                if len(fields) == 4:
                    v = [ffloat(i) for i in fields[1:]]
                else:
                    v = []
                eband.append(e)
                vband.append(v)
            ebands1.append(eband)
            mommat1.append(vband)
    kpoints = np.array(kpoints)
    ebands1 = np.array(ebands1)
    mommat1 = np.array(mommat1)
    # When several spin channels are present, the full list of k points is
    # redundant.
    kpoints = kpoints[:kpoints.shape[0] // nspin, :]

    # Lists of bands for different spin channels at the same k point are
    # concatenated.
    ebands = np.empty((nk, nspin, minband))
    for ispin in range(nspin):
        for ik in range(nk):
            ebands[ik, ispin, :] = ebands1[ispin * nk + ik][:minband]
    ebands = ebands.reshape((nk, nspin * minband))
    if mommat1.ndim == 3 and mommat1.shape[2] == 3:
        mommat = np.empty((nk, nspin, minband, 3))
        for ispin in range(nspin):
            for ik in range(nk):
                mommat[ik, ispin, :, :] = mommat1[ispin * nk + ik][:minband]
        mommat = mommat.reshape((nk, nspin * minband, 3))
    else:
        mommat = None

    # Convert everything to Ha and Ha / bohr
    efermi *= .5
    ebands *= .5
    if mommat is not None:
        mommat *= .5
    if nspin == 1:
        dosweight = 2.
    else:
        dosweight = 1.
    return efermi, dosweight, kpoints, ebands.T, mommat


def read_GENE_struct(filename):
    """Read in a GENE .structure file.

    Such files contain lattice vectors, elements and atomic positions.

    Args:
        filename: path to the .structure file
          NB: naming was previously .struct, and now .structure to distinguish 
          a complete file with all information needed for ASE
    Returns:
        An ASE Atoms object with crystal unit cell, 3D periodic, and appropriate atoms and positions
    """
    lines = open(filename, "r").readlines()
    # Line 1: title
    # Line 2-4: lattice vectors in bohr
    cell = np.empty((3, 3))
    for i in range(3):
        cell[i, :] = np.array([ffloat(j) for j in lines[i + 1].split()])
    # Line 5: number of atoms
    natom = int(lines[4].split()[0])
    # All remaining lines: element names and Cartesian coordinates in bohr
    cart = []
    symbolstring = ""
    for i in range(natom):
        fields = lines[5 + i].split()
        symbolstring += fields[0]
        cart.append([ffloat(j) for j in fields[1:4]])
    cart = np.array(cart)
    # Bundle everything in an Atoms object and return it
    return ase.Atoms(
        symbolstring,
        positions=cart / Angstrom,
        cell=cell / Angstrom,
        pbc=[1, 1, 1])


def W2Kmommat(filename, kpoints):
    """Read the contents of a Wien2k .mommat2 file.

    Args:
        filename: path to the .mommat2 file.
        kpoints: array with the k points of the same system.

    Returns:
        A 3-tuple. The first element is an array with the momentum matrix
        elements. The second and third are integers delimiting which bands
        have derivative information available.
    """
    nk = len(kpoints)

    fmat = open(filename, "r")
    matlines = fmat.readlines()
    fmat.close()
    il = 2
    mommat = []
    brk = []
    for ik in range(nk):
        nemin = int(matlines[il].split()[5])
        nemax = int(matlines[il].split()[6])
        il += 2
        brk += [[nemin, nemax]]
        ne = nemax - nemin + 1
        mmk = np.zeros((ne, 3), dtype=complex)
        ib1 = 0
        for ib in range(ne):
            line = matlines[il + ib1]
            for i in range(3):
                mmk[ib, i] = complex(
                    float(line[11 + 2 * i * 13:24 + 2 * i * 13]),
                    float(line[24 + 2 * i * 13:37 + 2 * i * 13]))
            ib1 += ne - ib
        mommat += [mmk]
        il += ib1 + 1
    nemax = np.min(brk, axis=0)[1]
    nemin = np.max(brk, axis=0)[0]
    mommat2 = []
    for ik in range(nk):
        mommat2 += [mommat[ik][nemin - brk[ik][0]:nemax - brk[ik][0] + 1]]
    return np.array(mommat2), nemin, nemax


def W2Kene(filename, conv):
    """Read the contents of a Wien2k .energy file.

    Args:
        filename: path to the .energy file.
        conv: result of passing the lattice vectors of the system through
            ase.io.wien2k_c2p.

    Returns:
        A 2-tuple. The first element is an (nkpoints, 3) array with the
        coordinates of the irreducible k points of the system. The second
        contains the energy bands in an array.
    """
    f = open(filename, "r", encoding="ascii")
    lines = f.readlines()
    f.close()
    linenumber = 0
    k = np.zeros(3)
    while True:
        try:
            ll = lines[linenumber]
            k[0], k[1], k[2], nband = float(ll[:19]), float(ll[19:38]), float(
                ll[38:57]), int(ll[73:79])
            break
        except:
            linenumber += 1

    minband = np.infty
    ebands1 = []
    kpoints = []
    while True:
        try:
            ll = lines[linenumber]
            k[0], k[1], k[2], nband = float(ll[:19]), float(ll[19:38]), float(
                ll[38:57]), int(ll[73:79])
            nband = int(nband)
            if nband < minband:
                minband = nband
            linenumber += 1
            eband = []
            for i in range(nband):
                e = ffloat(lines[linenumber].split()[1])
                eband += [e]
                linenumber += 1
            ebands1 += [np.array(eband)]
            kpoints += [k.copy()]
        except:
            break
    kpoints = kpoints @ conv
    ebands = np.zeros((len(kpoints), minband))
    for i in range(len(kpoints)):
        ebands[i] = ebands1[i][:minband]
    return kpoints, ebands.T * .5  # To Ha


def W2Kfermi(filename):
    """Read the value of the Fermi level from a Wien2k .scf file.

    Args:
        filename: path to the .scf file.

    Returns:
        The value of the Fermi level.
    """
    with open(filename) as f:
        for l in f:
            if l.startswith(":FER"):
                return .5 * float(l[38:53])
    return None


# Small utility functions for parsing VASP's xml files.
def _parse_vasp_array(line, cast=float):
    """Transform a line of text with fields separated by spaces into a list
    of homogeneous objects.
    """
    return [cast(i) for i in line.split()]


# Specialized parsers for each piece of information
def _parse_vasp_name(xml_tree):
    """Extract the system name from an XML ElementTree representing a
    vasprun.xml file.
    """
    xml_general = xml_tree.find("./parameters/separator[@name=\"general\"]")
    return xml_general.find("./i[@name=\"SYSTEM\"]").text


def _parse_vasp_fermi(xml_tree):
    """Extract the Fermi level from an XML ElementTree representing a
    vasprun.xml file.
    """
    xml_dos = xml_tree.find("./calculation/dos")
    return float(xml_dos.find("./i[@name=\"efermi\"]").text)


def _parse_vasp_nelect(xml_tree):
    """Extract NELECT from an XML ElementTree representing a vasprun.xml file.
    """
    xml_electronic = xml_tree.find("./parameters").find(
        "./separator[@name=\"electronic\"]")
    return float(xml_electronic.find("./i[@name=\"NELECT\"]").text)


def _parse_vasp_structure(xml_tree):
    """Extract the structural information from an XML ElementTree representing a
    vasprun.xml file, and return it as an ASE atoms object.
    """
    # Read in the lattice vectors, reduced positions and chemical elements.
    lattvec = np.empty((3, 3))
    positions = []
    elements = []
    xml_basis = xml_tree.find(
        "./calculation/structure/crystal/varray[@name=\"basis\"]")
    for i, line in enumerate(xml_basis):
        lattvec[i, :] = _parse_vasp_array(line.text)
    xml_positions = xml_tree.find(
        "./calculation/structure/varray[@name=\"positions\"]")
    for line in xml_positions:
        positions.append(_parse_vasp_array(line.text))
    positions = np.array(positions)
    xml_elements = xml_tree.findall("./atominfo/array[@name=\"atoms\"]/set/rc")
    for element in xml_elements:
        elements.append(element.find("./c").text)
    # Build and return an ase atoms object
    cartesian = positions @ lattvec
    atoms = ase.Atoms(
        elements, positions=cartesian, cell=lattvec, pbc=[1, 1, 1])
    return atoms


def _parse_vasp_ikpoints(xml_tree):
    """Extract a list of irreducible k points from an XML ElementTree
    representing a vasprun.xml file, and return it as an (nkpoints, 3) numpy
    array.
    """
    xml_ikpoints = xml_tree.find("./kpoints/varray[@name=\"kpointlist\"]")
    ikpoints = []
    for p in xml_ikpoints:
        ikpoints.append(_parse_vasp_array(p.text))
    return np.array(ikpoints)


def _parse_vasp_eigenvalues(xml_tree):
    """Extract a list of eigenvalues at each irreducible k point from an XML
    ElementTree representing a vasprun.xml file, and return it as an
    (nspin, nkpoints, nbands) numpy array.
    """
    xml_eigenvalues = xml_tree.find("./calculation/eigenvalues/array/set")
    data = []
    for spin in xml_eigenvalues:
        dspin = []
        for point in spin:
            dpoint = []
            for band in point:
                # Each row contains four fields: energy and occupancy
                dpoint.append(_parse_vasp_array(band.text)[0])
            dspin.append(dpoint)
        data.append(dspin)
    return np.array(data)


def _parse_vasp_velocities(xml_tree):
    """Extract a list of k points and a list of group velocities at those
    points from an XML ElementTree representing a vasprun.xml file. Return
    them as a tuple of two numpy arrays, with shapes (nkpoints, 3) and
    (nspin, nkpoints, nbands, 3). Note that these k points are not symmetry
    reduced. Return None if the tree does not contain information about group
    velocities.
    """
    xml_vels = xml_tree.find("*/electronvelocities")
    if xml_vels is None:
        return None
    xml_kpoints = xml_vels.find("./kpoints/varray")
    kpoints = []
    for p in xml_kpoints:
        kpoints.append(_parse_vasp_array(p.text))
    # The three axes of the data are ordered like in (spin, point, band),
    # from most to least significant.
    xml_data = xml_vels.find("./eigenvalues/array/set")
    data = []
    for spin in xml_data:
        dspin = []
        for point in spin:
            dpoint = []
            for band in point:
                # Each row contains four fields: energy, vx, vy, vz
                dpoint.append(_parse_vasp_array(band.text)[1:])
            dspin.append(dpoint)
        data.append(dspin)
    return (np.array(kpoints), np.array(data))


def parse_vasprunxml(filename):
    """Parse a vasprun.xml file to extract all structural information, the
    coordinates of the q points used in the calculation, the energy eigenvalues
    at those points, and their gradients (group velocities) if available.

    For vasp to save the group velocities, both LOPTICS and LVEL must have been
    set to .TRUE. in the INCAR file.

    Return the results in a dictionary with the following keys:
    name = the system name (often listed as "unknown system" in vasprun.xml)
    atoms: an ASE atoms object representing the structure
    nelect: number of valence electrons in the simulation cell
    kpoints: numpy array with shape (nk, 3) containing direct coordinates of an
             irreducible set of k points.
    fermi: the Fermi level in eV
    E: numpy array with shape (nspin, npoints, nbands) containing the energies
       of the electronic eigenstates at each point
    v: numpy array with shape (nspin, npoints, nbands, 3) containing the
       gradient of E with respect to k at each point. This key is only present
       if the vasprun.xml file contains group velocities.
    """
    # Open the file and parse it with the builtin parser
    xml_tree = et.parse(filename)
    # Obtain each piece of information from the corresponding function.
    nruter = dict(
        fermi=_parse_vasp_fermi(xml_tree),
        name=_parse_vasp_name(xml_tree),
        atoms=_parse_vasp_structure(xml_tree),
        nelect=_parse_vasp_nelect(xml_tree),
        kpoints=_parse_vasp_ikpoints(xml_tree),
        E=_parse_vasp_eigenvalues(xml_tree))
    # Group velocities require special care. They may not be present, and
    # even if they are we need to extract the subset of data corresponding
    # to the irreducible k points.
    res = _parse_vasp_velocities(xml_tree)
    if res is not None:
        k, v = res
        # Obtain the index of those k points from the full array that are also
        # in the array of irreducible k points.
        reduk = nruter["kpoints"] % 1.
        fullk = k % 1.
        distances = sp.spatial.distance.cdist(reduk, fullk)
        indices = []
        for i in range(reduk.shape[0]):
            pos = distances[i, :].argmin()
            if distances[i, pos] > 1e-6:
                raise ValueError("inconsistent sets of k points")
            indices.append(pos)
        if len(set(indices)) != len(indices):
            raise ValueError("inconsistent sets of k points")
        nruter["v"] = v[:, indices, :, :]
    return nruter


def save_trace(filename, data, Tr, mur, N, sdos, cv, cond, seebeck, kappa,
               hall):
    """Create a .trace file in a format similar to the original output of
    BoltzTraP.
    """
    nformulas = data.get_formula_count()
    headerfmt = "#{:>9s} {:>9s}" + " ".join(8 * ["{:>25s}"])
    header = [""] * 10
    header[0] = "Ef[Ry]"
    header[1] = "T[K]"
    header[2] = "N[e/uc]"
    header[3] = "DOS(ef)[1/(Ha*uc)]"
    header[4] = "S[V/K]"
    header[5] = "sigma/tau0[1/(ohm*m*s)]"
    header[6] = "RH[m**3/C]"
    header[7] = "kappae/tau0[W/(m*K*s)]"
    header[8] = "cv[J/(mol*K)]"
    header[9] = "chi[m**3/mol]"
    rowfmt = "{:>10g} {:>9g}" + " ".join(8 * ["{:>25g}"])
    with open(filename, "w") as f:
        print(headerfmt.format(*header).strip(), file=f)
        for imu, mu in enumerate(mur):
            for iT, T in enumerate(Tr):
                # cv is expressed as a quantity per unit cell. To create the
                # output, we reexpress it as a quantity per mole of unit
                # formula.
                ocv = cv[iT, imu] * AVOGADRO / nformulas
                # The equivalent to the trace of the Hall tensor is an average
                # over the even permutations of [0, 1, 2].
                ohall = (hall[iT, imu, 0, 1, 2] + hall[iT, imu, 2, 0, 1] +
                         hall[iT, imu, 1, 2, 0]) / 3.
                # Our estimate of the susceptibility comes directly from the
                # smoothed DOS.
                magsus = sdos[iT, imu] * MU0 * MUB**2 * AVOGADRO / (
                    data.atoms.get_volume() * Meter**3)
                print(
                    rowfmt.format(mu * 2., T, N[iT, imu] - data.nelect,
                                  sdos[iT, imu], seebeck[iT, imu].trace() / 3.,
                                  cond[iT, imu].trace() / 3., ohall,
                                  kappa[iT, imu].trace() / 3., ocv, magsus),
                    file=f)


def save_condtens(filename, Tr, mur, N, cond, seebeck, kappa):
    """Create a .condtens file in a format similar to the original output of
    BoltzTraP.
    """
    headerfmt = "#{:>9s} {:>9s}" + " ".join(28 * ["{:>25s}"])
    header = [""] * 30
    header[0] = "Ef[Ry]"
    header[1] = "T[K]"
    header[2] = "N[e/uc]"
    header[3] = "sigma/tau0[1/(ohm*m*s)]"
    header[12] = "S[V/K]"
    header[21] = "kappae/tau0[W/(m*K*s)]"
    rowfmt = "{:>10g} {:>9g}" + " ".join(28 * ["{:>25g}"])
    with open(filename, "w") as f:
        print(headerfmt.format(*header).strip(), file=f)
        for imu, mu in enumerate(mur):
            for iT, T in enumerate(Tr):
                print(
                    rowfmt.format(
                        mu * 2.,
                        T,
                        N[iT, imu],
                        *cond[iT, imu, :, :].ravel(order="F"),
                        *seebeck[iT, imu, :, :].ravel(order="F"),
                        *kappa[iT, imu, :, :].ravel(order="F")),
                    file=f)


def save_halltens(filename, Tr, mur, N, hall):
    """Create a .halltensfile in a format similar to the original output of
    BoltzTraP.
    """
    headerfmt = "#{:>9s} {:>9s}" + " ".join(28 * ["{:>25s}"])
    header = [""] * 30
    header[0] = "Ef[Ry]"
    header[1] = "T[K]"
    header[2] = "N[e/uc]"
    header[3] = "RH[m**3/C]"
    rowfmt = "{:>10g} {:>9g}" + " ".join(28 * ["{:>25g}"])
    with open(filename, "w") as f:
        print(headerfmt.format(*header).strip(), file=f)
        for imu, mu in enumerate(mur):
            for iT, T in enumerate(Tr):
                print(
                    rowfmt.format(
                        mu * 2.,
                        T,
                        N[iT, imu],
                        *hall[iT, imu, :, :, :].ravel(order="F")),
                    file=f)
