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

import io
import os
import os.path
import sys
import glob
import logging
import math
import collections
import functools

import numpy as np
import ase.io
import ase.io.wien2k

import BoltzTraP2.io
import BoltzTraP2.misc
import BoltzTraP2.sphere
from BoltzTraP2.units import *


class LoaderError(Exception):
    """Generic exception class for problems in loaders."""
    pass


# List of format names and loaders that will be tried for each directory
# A loader can be any callable that takes a mandatory argument (usually a
# directory) and returns an object with the following mandatory attributes:
#
# sysname: an arbitrary string
# atoms: an ASE atoms object representing the system
# dosweight: maximum occupancy of each state (normally 1 for spin-polarized
#            calculations and 2 otherwise)
# kpoints: irreducible set of k points in reduced coordinates, as a
#          (nkpoints, 3) array
# fermi: Fermi level from the original calculation, in Ha
# ebands: energy bands in Ha, as a (nbands, nkpoints) array
#
# And some or all of the following optinal attributes:
#
# self.nelect: number of valence electrons in the calculation
# self.mommat: momentum matrix elements in in Ha / bohr, as a
#              (nkpoints, nbands, 3) array
#
# The last loader to be added has the highest priority.
# The loader must raise a LoaderException if the argument cannot be processed.
# Note that the argument can be of any type. If the loader expects a string,
# it must explicitly check that fact.
# Implementing a new loader and adding it to this list is the
# way to add support for a new format in BoltzTraP2.
# To add a loader to the tuple, use BolzTraP2.dft.register_loader(name, loader)
# Example: register_loader("VASP", VASPLoader)
loaders = []


def register_loader(name, loader):
    """Add an element to the list of registered loaders."""
    loaders.append((name, loader))


class VASPLoader:
    """Loader for VASP calculations."""

    def __init__(self, directory):
        if not isinstance(directory, str):
            raise LoaderError("this loader only works with directories")
        with BoltzTraP2.misc.dir_context(directory):
            vasprunxml = "vasprun.xml"
            if not os.path.isfile(vasprunxml):
                raise LoaderError("vasprun.xml not found")
            r = BoltzTraP2.io.parse_vasprunxml(vasprunxml)
            self.sysname = r["name"]
            self.atoms = r["atoms"]
            BoltzTraP2.misc.info("lattice:", self.atoms.get_cell().T)
            self.kpoints = r["kpoints"]
            # Some unit conversions are necessary at this point.
            self.fermi = r["fermi"] * eV
            # If the calculation is spin-polarized, concatenate all the bands
            self.ebands = r["E"].transpose((0, 2, 1)) * eV
            if self.ebands.shape[0] > 1:
                self.dosweight = 1.0
            else:
                self.dosweight = 2.0
            self.ebands = self.ebands.reshape(-1, self.ebands.shape[-1])
            self.nelect = r["nelect"]
            if "v" in r:
                self.mommat = r["v"].transpose((0, 2, 1, 3)) * eV / Angstrom
                self.mommat = self.mommat.reshape(
                    tuple([-1] + list(self.mommat.shape[2:])))
                self.mommat = self.mommat.transpose((1, 0, 2))


register_loader("VASP", VASPLoader)


class GenericWien2kLoader:
    """Generic Loader for Wien2k calculations, not intended for direct use.

    This class gives total control over the file names to be loaded.
    """

    def __init__(self, w2kname, dosweight, scffn, structfn, energyfn,
                 derfn=""):
        BoltzTraP2.misc.info("Wien2k system name:", w2kname)
        self.sysname = w2kname
        if not os.path.isfile(scffn):
            raise LoaderError(".scf file not found")
        self.dosweight = dosweight
        if not os.path.isfile(energyfn):
            raise LoaderError("energy(so) file not found")
        BoltzTraP2.misc.info("matching energy and scf files were found")
        self.atoms = ase.io.read(structfn)
        dum = ase.io.wien2k.read_struct(structfn, ase=False)
        latt = dum[1]
        if latt == "R":
            latt = "P"
        conv = ase.io.wien2k.c2p(latt)
        BoltzTraP2.misc.info("lattice:", self.atoms.get_cell().T)
        BoltzTraP2.misc.info("conv:", conv)
        BoltzTraP2.misc.info("conv:", latt)
        self.kpoints, self.ebands = BoltzTraP2.io.W2Kene(energyfn, conv)
        self.fermi = BoltzTraP2.io.W2Kfermi(scffn)
        if os.path.isfile(derfn):
            BoltzTraP2.misc.info("a matching .mommat2 file was found")
            self.mommat, nemin, nemax = BoltzTraP2.io.W2Kmommat(derfn,
                                                                self.kpoints)
            self.ebands = self.ebands[nemin - 1:nemax]


class Wien2kLoader(GenericWien2kLoader):
    """Loader for Wien2k calculations following the usual naming convention."""

    def __init__(self, directory):
        """Guess the name of the system and call the base cosntructor to do
        the real work.
        """
        if not isinstance(directory, str):
            raise LoaderError("this loader only works with directories")
        w2kname = _get_W2Ksystemname(directory)
        if w2kname is None:
            raise LoaderError("cannot determine a Wien2k system name")
        scffn = w2kname + ".scf"
        structfn = w2kname + ".struct"
        energyfn = w2kname + ".energyso"
        with BoltzTraP2.misc.dir_context(directory):
            if not os.path.isfile(energyfn):
                energyfn = os.path.splitext(energyfn)[0] + ".energy"
        if os.path.splitext(energyfn)[1] == ".energy":
            dosweight = 2.
        else:
            dosweight = 1.
        derfn = w2kname + ".mommat2"
        GenericWien2kLoader.__init__(self, w2kname, dosweight, *(
            os.path.join(directory, i)
            for i in (scffn, structfn, energyfn, derfn)))


register_loader("Wien2k", Wien2kLoader)


class GENELoader:
    """Loader for data in BoltzTraP's generic format."""

    def __init__(self, directory):
        if not isinstance(directory, str):
            raise LoaderError("this loader only works with directories")
        genename = _get_GENEsystemname(directory)
        if genename is None:
            raise LoaderError("cannot determine a GENE system name")
        structfn = genename + ".structure"
        energyfn = genename + ".energy"
        with BoltzTraP2.misc.dir_context(directory):
            BoltzTraP2.misc.info("GENE system name:", genename)
            if not os.path.isfile(energyfn):
                raise ValueError("energy file not found")
            self.atoms = BoltzTraP2.io.read_GENE_struct(structfn)
            BoltzTraP2.misc.info("lattice:", self.atoms.get_cell().T)
            (self.fermi, self.dosweight, self.kpoints, self.ebands,
             mommat) = BoltzTraP2.io.read_GENE_eneandmat(energyfn)
            if mommat is not None:
                self.mommat = mommat
        self.sysname = genename


register_loader("GENE", GENELoader)


def _get_GENEsystemname(dirname):
    """Try to guess the GENE system name corresponding to a directory."""
    with BoltzTraP2.misc.dir_context(dirname):
        filenames = sorted(
            [i for i in glob.glob("*.structure") if os.path.isfile(i)])
        if not filenames:
            return None
        if len(filenames) > 1:
            logging.warning(
                "there is more than one .structure file in the directory"
                " - using the first one")
    return os.path.splitext(os.path.basename(filenames[0]))[0]


def _get_W2Ksystemname(dirname):
    """Try to guess the Wien2k system name corresponding to a directory."""
    with BoltzTraP2.misc.dir_context(dirname):
        filenames = sorted(
            [i for i in glob.glob("*.struct") if os.path.isfile(i)])
        if not filenames:
            return None
        if len(filenames) > 1:
            logging.warning(
                "there is more than one .struct file in the directory "
                "- using the first one")
    return os.path.splitext(os.path.basename(filenames[0]))[0]


class DFTData:
    """Objects of this class hold structural and dynamical information from DFT
    results in any supported format.
    """

    def __init__(self, directory, derivatives=False, *args, **kwargs):
        """Create a DFTData object."""
        for label, loader in loaders[::-1]:
            BoltzTraP2.misc.info("looking for a {} calculation".format(label))
            try:
                loaded = loader(directory, *args, **kwargs)
            except LoaderError as e:
                BoltzTraP2.misc.info("error in {} loader: {}".format(label, e))
                continue
            self.source = label
            break
        else:
            raise ValueError(
                "no calculation found in directory {}".format(directory))
        BoltzTraP2.misc.info(
            "successfully loaded a {} calculation".format(self.source))
        # Try to copy all relevant attributes from the loader
        if derivatives:
            try:
                self.mommat = loaded.mommat
            except AttributeError:
                raise ValueError(
                    "no derivative information found in directory {}".format(
                        directory))
        else:
            try:
                loaded.mommat
            except AttributeError:
                pass
            else:
                BoltzTraP2.misc.info(
                    "derivative information will be discarded")
            self.mommat = None
        try:
            self.sysname = loaded.sysname
            self.atoms = loaded.atoms
            self.dosweight = loaded.dosweight
            self.kpoints = loaded.kpoints
            self.fermi = loaded.fermi
            self.ebands = loaded.ebands
        except AttributeError:
            raise ValueError(
                "some essential piece of information was not loaded")
        BoltzTraP2.misc.info("Fermi energy:", self.fermi)
        # If the number of valence electrons has not been set yet, compute
        # it from the bands
        try:
            self.nelect = loaded.nelect
        except AttributeError:
            degeneracies = BoltzTraP2.sphere.calc_reciprocal_degeneracies(
                self.atoms, self.kpoints)
            weights = degeneracies.astype(np.float64) / degeneracies.sum()
            occupancy = (loaded.ebands < loaded.fermi).astype(np.int)
            self.nelect = round(self.dosweight * (occupancy * weights).sum())

    def bandana(self, emin=-np.inf, emax=np.inf):
        bandmin = np.min(self.ebands, axis=1)
        bandmax = np.max(self.ebands, axis=1)
        II = np.nonzero(bandmin < emax)
        nemax = II[0][-1]
        II = np.nonzero(bandmax > emin)
        nemin = II[0][0]
        BoltzTraP2.misc.info("BANDANA output")
        for iband in range(len(self.ebands)):
            BoltzTraP2.misc.info(iband, bandmin[iband], bandmax[iband], (
                (bandmin[iband] < emax) & (bandmax[iband] > emin)))
        self.ebands = self.ebands[nemin:nemax]
        if self.mommat is not None:
            self.mommat = self.mommat[:, nemin:nemax, :]
        # Removing bands may change the number of valence electrons
        self.nelect -= self.dosweight * nemin
        return nemin, nemax

    def get_lattvec(self):
        try:
            self.lattvec
        except AttributeError:
            self.lattvec = self.atoms.get_cell().T * Angstrom
        return self.lattvec

    def get_volume(self):
        try:
            self.UCVol
        except AttributeError:
            lattvec = self.get_lattvec()
            self.UCvol = np.abs(np.linalg.det(lattvec))
        return self.UCvol

    def get_formula_count(self):
        """Return the number of irreducible formulas in the unit cell.

        Useful for computing molar quantities.
        """
        counts = collections.Counter(self.atoms.get_chemical_symbols())
        return functools.reduce(math.gcd, counts.values())
