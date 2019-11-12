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

import ast
import sys
import os
import os.path
import glob
import itertools
import logging
import argparse
import contextlib

import numpy as np
import scipy as sp
import scipy.linalg as la
import scipy.spatial
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import ase.dft.kpoints as asekp

try:
    import colorama
    colorama.init()
    terminal_colors = True
except ModuleNotFoundError:
    terminal_colors = False

import BoltzTraP2
import BoltzTraP2.dft
import BoltzTraP2.fite
import BoltzTraP2.sphere
import BoltzTraP2.bandlib
import BoltzTraP2.serialization
import BoltzTraP2.fermisurface
from BoltzTraP2.misc import info
from BoltzTraP2.misc import warning
from BoltzTraP2.misc import lexit
from BoltzTraP2.misc import dir_context
from BoltzTraP2.misc import TimerContext
from BoltzTraP2.version import *
from BoltzTraP2.units import *

PSEUDO_BLACK = "#333333"


def set_logging_level(verbose):
    """Set the logging level based on the desired verbosity"""
    logger = logging.getLogger()
    if verbose == 0:
        logger.setLevel(logging.ERROR)
    elif verbose == 1:
        logger.setLevel(logging.WARNING)
    elif verbose == 2:
        logger.setLevel(logging.INFO)
    elif verbose >= 3:
        logger.setLevel(logging.DEBUG)
    else:
        raise ValueError("the argument must be a positive integer")


def positive_integer(arg):
    """Cast arg to a positive integer or raise an exception."""

    def error():
        raise Vargparse.ArgumentTypeError(
            "{} is not a positive integer".format(arg))

    try:
        out = int(arg)
    except ValueError:
        error()
    if out < 1:
        error()
    return out


def positive_float(arg):
    """Cast arg to a positive floating-point number or raise an exception."""

    def error():
        raise argparse.ArgumentTypeError(
            "{} is not a positive number".format(arg))

    try:
        out = float(arg)
    except ValueError:
        error()
    if out < 0.:
        error()
    return out


def temperature_argument(arg):
    """Cast arg to a NumPy array of temperatures.

    Three formats are accepted:
    1. A single floating-point number
    2. A list of floating-point values
    3. A 3-tuple (Tmin, Tmax, nT) used to build a range
    """
    try:
        value = ast.literal_eval(arg)
    except:
        raise argparse.ArgumentTypeError(
            "'{}' cannot be parsed as a Python literal".format(arg))
    if isinstance(value, (int, float)):
        return np.array([float(value)])
    if isinstance(value, list):
        try:
            nruter = np.array(sorted(value), dtype=np.float64)
        except ValueError:
            raise argparse.ArgumentTypeError(
                "some temperatures cannot be cast to floats")
        return nruter
    if isinstance(value, tuple):
        if len(value) != 3:
            raise argparse.ArgumentTypeError(
                "a temperature range must be a 3-tuple")
        try:
            Tmin = float(value[0])
        except ValueError:
            raise argparse.ArgumentTypeError("Tmin cannot be cast to a float")
        try:
            Tmax = float(value[1])
        except ValueError:
            raise argparse.ArgumentTypeError("Tmax cannot be cast to a float")
        if Tmax <= Tmin:
            raise argparse.ArgumentTypeError("Tmax must be higher than Tmin")
        nT = value[2]
        if not isinstance(nT, int) or nT <= 0:
            raise argparse.ArgumentTypeError("nT is not a positive integer")
        return np.linspace(Tmin, Tmax, num=nT, endpoint=True)
    raise argparse.ArgumentTypeError("invalid temperature specification")


def components_argument(arg):
    """Parse an argument into a set of tensor components.

    The input must be a list in Python syntax, and each of its elements a
    string of length 2 or 3 made up of the letters x, y and z, or the special
    value "scalar" meaning 1/3 of the trace.
    """
    axes = "xyz"
    try:
        value = ast.literal_eval(arg)
    except:
        raise argparse.ArgumentTypeError(
            "'{}' cannot be parsed as a Python literal".format(arg))
    if not (isinstance(value, list) and all(isinstance(i, str)
                                            for i in value)):
        raise argparse.ArgumentTypeError(
            "'{}' cannot be parsed as a Python list of strings".format(arg))
    # Remove duplicates before parsing
    value = set(value)
    nruter = []
    for i in value:
        if i == "scalar":
            nruter.append(None)
        elif len(i) in (2, 3):
            component = []
            for j in i:
                try:
                    component.append(axes.index(j))
                except ValueError:
                    raise argparse.ArgumentTypeError(
                        "character '{}' is not valid in a tensor component"
                        " specification".format(j))
            nruter.append(component)
        else:
            raise argparse.ArgumentTypeError(
                "tensor component specifications can only contain strings of"
                " length 2 or 3, or the special value 'scalar'")
    return nruter


def is_writable(filename):
    """Check if a file exists and is writeable, or if it could be created."""
    if os.path.isfile(filename) and os.access(filename, os.W_OK):
        return True
    dirname = os.path.abspath(os.path.dirname(filename))
    return os.path.isdir(dirname) and os.access(dirname, os.W_OK)


def parse_interpolate(args):
    """Parse command-line arguments for the interpolate subcommand."""
    # Check if the output file is writable at this point. Note that this
    # does not guarantee that it remains so. It is just a "courtesy" check
    # to avoid a lengthy calculation ending in failure to save.
    args.directory = os.path.abspath(args.directory)
    args.output = os.path.abspath(args.output)
    if not os.path.isdir(args.directory):
        lexit("the specified directory does not exist")
    if not is_writable(args.output):
        lexit("the output file is not writable")
    # Try and read the input.
    data = BoltzTraP2.dft.DFTData(args.directory, args.derivatives)
    # Perform some sanity checks on the energy window specification
    if args.emin >= args.emax:
        lexit("zero-width energy window")
    emin = args.emin + (0. if args.absolute else data.fermi)
    emax = args.emax + (0. if args.absolute else data.fermi)
    if emin >= data.fermi or emax <= data.fermi:
        lexit("the energy window must bracked the Fermi level")
    info("Nvalence (before BANDANA) =", data.nelect)
    # Drop bands outside the chosen energy range
    nemin = data.bandana(emin, emax)[0]
    info("Nvalence (after BANDANA) =", data.nelect)
    # Refuse to interpolate to fewer k points than there are in the input
    nkinput = data.kpoints.shape[0]
    logging.info(
        "there are {} irreducible k points in the input".format(nkinput))
    if args.multiplier is not None:
        nktarget = args.multiplier * nkinput
    else:
        nktarget = args.kpoints
    logging.info(
        "{} irreducible k points have been requested".format(nktarget))
    equivalences = BoltzTraP2.sphere.get_equivalences(data.atoms, nktarget)
    nkoutput = len(equivalences)
    logging.info(
        "{} irreducible k points have been generated".format(nkoutput))
    if nkoutput < nkinput:
        lexit("refusing to interpolate to a sparser grid")
    # Perform the interpolation
    with TimerContext() as timer:
        coeffs = BoltzTraP2.fite.fitde3D(data, equivalences, args.nworkers)
        deltat = timer.get_deltat()
        info("the interpolation took {:.3g} s".format(deltat))
    # Gather the metadata
    metadata = BoltzTraP2.serialization.gen_bt2_metadata(data,
                                                         args.derivatives)
    # Save the result
    info("about to save the results to", args.output)
    with TimerContext() as timer:
        BoltzTraP2.serialization.save_calculation(
            args.output, data, equivalences, coeffs, metadata)
        deltat = timer.get_deltat()
        info("saving the results took {:.3g} s".format(deltat))


def parse_integrate(args):
    """Parse command-line arguments for the integrate subcommand."""
    # If the number of bins is not set by hand, let the code handle it
    # automatically.
    try:
        args.bins
    except AttributeError:
        args.bins = None
    # Try and load the data from the interpolation step
    data, equivalences, coeffs, metadata = (
        BoltzTraP2.serialization.load_calculation(args.bt2_file))
    lattvec = data.get_lattvec()
    info("sucessfully loaded " + args.bt2_file)
    # Rebuild the bands
    with TimerContext() as timer:
        eband, vvband, cband = BoltzTraP2.fite.getBTPbands(
            equivalences, coeffs, lattvec, True, args.nworkers)
        deltat = timer.get_deltat()
        info("rebuilding the bands took {:.3g} s".format(deltat))

    # Compute the DOS histogram
    with TimerContext() as timer:
        epsilon, dos, vvdos, cdos = BoltzTraP2.bandlib.BTPDOS(
            eband,
            vvband,
            cband,
            npts=args.bins,
            scattering_model=args.scattering_model)
        deltat = timer.get_deltat()
        info("computing the DOS took {:.3g} s".format(deltat))
    info("Number of DOS bins:", epsilon.size)

    # Refine the estimate of the intrinsic chemical potential at
    # each temperature.
    Tr = args.temperature
    mu0 = np.empty_like(Tr)
    for iT, T in enumerate(Tr):
        mu0[iT] = BoltzTraP2.bandlib.refine_mu0(epsilon, dos, data.nelect, T,
                                                data.dosweight)
    # And at 0 K
    fermi = BoltzTraP2.bandlib.refine_mu0(epsilon, dos, data.nelect, 0.,
                                          data.dosweight)

    # Compute the moments of the FD distribution at each point.
    # Determine the chemical potentials to explore by taking the values of the
    # energy bins and discarding those that are too close to the edge to give
    # meaningful results.
    # 9 * kB * T gives a reduction in fFD of about 1e-4, a relatively
    # safe margin (although still far from the hard cutoff in bandlib).
    margin = 9. * BOLTZMANN * Tr.max()
    mur_indices = np.logical_and(epsilon > epsilon.min() + margin,
                                 epsilon < epsilon.max() - margin)
    mur = epsilon[mur_indices]
    if mur.size == 0:
        lexit("the energy window is too narrow")

    # Get the smooth DOS at each temperature
    with TimerContext() as timer:
        sdos = np.empty((Tr.size, epsilon.size))
        for iT, T in enumerate(Tr):
            sdos[iT, :] = BoltzTraP2.bandlib.smoothen_DOS(epsilon, dos, T)
        sdos = sdos[:, mur_indices]
        deltat = timer.get_deltat()
        info("smoothing the DOS took {:.3g} s".format(deltat))

    info("Temperatures:", Tr)
    info("Fermi level from DFT:", data.fermi)
    info("Refined Fermi level:", fermi)
    info("Intrinsic chemical potential for each T:", mu0)
    info("Chemical potentials:", mur)
    with TimerContext() as timer:
        cv = BoltzTraP2.bandlib.calc_cv(
            epsilon, dos, mur, Tr, dosweight=data.dosweight)
        N, L0, L1, L2, Lm11 = BoltzTraP2.bandlib.fermiintegrals(
            epsilon,
            dos,
            vvdos,
            mur=mur,
            Tr=Tr,
            dosweight=data.dosweight,
            cdos=cdos)
        deltat = timer.get_deltat()
        info("computing the FD moments took {:.3g} s".format(deltat))

    # Rescale and combine the moments to get the Onsager transport coefficients
    vuc = data.atoms.get_volume() * Angstrom**3
    L11, seebeck, kappa, hall = BoltzTraP2.bandlib.calc_Onsager_coefficients(
        L0, L1, L2, mur, Tr, vuc, Lm11)

    # Save the results to BoltzTraP-style files
    basefn = os.path.splitext(args.bt2_file)[0]
    tracefn = basefn + ".trace"
    info("trace output file:", tracefn)
    BoltzTraP2.io.save_trace(tracefn, data, Tr, mur, N, sdos, cv, L11, seebeck,
                             kappa, hall)
    condtensfn = basefn + ".condtens"
    info("conductivity/seebeck output file:", condtensfn)
    BoltzTraP2.io.save_condtens(condtensfn, Tr, mur, N, L11, seebeck, kappa)
    halltensfn = basefn + ".halltens"
    info("Hall coefficient output file:", halltensfn)
    BoltzTraP2.io.save_halltens(halltensfn, Tr, mur, N, hall)

    # Save the results to a single compressed JSON file
    metadata2 = BoltzTraP2.serialization.gen_btj_metadata(
        metadata, args.scattering_model)
    jsonfn = basefn + ".btj"
    info("JSON output file:", jsonfn)
    BoltzTraP2.serialization.save_results(jsonfn, data, fermi, Tr, mu0, mur, N,
                                          sdos, cv, L11, seebeck, kappa, hall,
                                          metadata2)


def parse_plotbands(args):
    """Parse command-line arguments for the plotbands subcommand."""
    # Try and load the data from the interpolation step
    data, equivalences, coeffs, metadata = (
        BoltzTraP2.serialization.load_calculation(args.bt2_file))
    lattvec = data.get_lattvec()
    info("sucessfully loaded " + args.bt2_file)
    # The second position alargument is first interpreted as a Python literal,
    # and after parsing it is cast to a NumPy array, which must have the right
    # dimensions. The special value None directs the parser to split the path
    # in several parts.
    try:
        kpaths = ast.literal_eval(args.kpath)
    except ValueError:
        lexit("'{}' cannot be parsed as a Python literal".format(kpaths))
    if not isinstance(kpaths, list):
        lexit("'{}' cannot be parsed as a Python list".format(kpaths))
    kpaths = [
        list(group)
        for k, group in itertools.groupby(
            kpaths, key=lambda x: x is not None) if k
    ]
    try:
        kpaths = [np.array(i, dtype=np.float64) for i in kpaths]
        for i in kpaths:
            if i.shape[0] < 2 or i.shape[1] != 3:
                raise ValueError
    except ValueError:
        lexit("the path cannot be interpreted as a set of N x 3"
              " arrays (with N >= 2")

    plt.figure()
    ax = plt.gca()
    ticks = []
    dividers = []
    offset = 0.
    for ikpath, kpath in enumerate(kpaths):
        ax.set_prop_cycle(
            color=matplotlib.rcParams["axes.prop_cycle"].by_key()["color"])
        info("k path #{}".format(i + 1))
        # Generate the explicit point list.
        kp, dkp, dcl = asekp.bandpath(kpath, data.atoms.cell, args.nkpoints)
        dkp += offset
        dcl += offset
        # Compute the band energies.
        with TimerContext() as timer:
            egrid = BoltzTraP2.fite.getBands(kp, equivalences,
                                             data.get_lattvec(), coeffs)[0]
            deltat = timer.get_deltat()
            info("rebuilding the bands took {:.3g} s".format(deltat))
        egrid -= data.fermi
        # Create the plot
        nbands = egrid.shape[0]
        for i in range(nbands):
            plt.plot(dkp, egrid[i, :], lw=2.)
        ticks += dcl.tolist()
        dividers += [dcl[0], dcl[-1]]
        offset = dkp[-1]
    ax.set_xticks(ticks)
    ax.set_xticklabels([])
    for d in ticks:
        plt.axvline(x=d, color=PSEUDO_BLACK, ls="--", lw=.5)
    for d in dividers:
        plt.axvline(x=d, color=PSEUDO_BLACK, ls="-", lw=2.)
    plt.axhline(y=0., color=PSEUDO_BLACK, lw=1.)
    plt.ylabel(r"$\varepsilon - \varepsilon_F\;\left[\mathrm{Ha}\right]$")
    plt.tight_layout()
    plt.show()


def parse_plot(args):
    """Parse command-line arguments for the plotbands subcommand."""
    # Try and load the data from the integration step
    (data, fermi, Tr, mu0, mur, N, sdos, cv, cond, seebeck, kappa, hall,
     metadata) = BoltzTraP2.serialization.load_results(args.btj_file)
    info("sucessfully loaded " + args.btj_file)
    nformulas = data.get_formula_count()
    # Perform some sanity checks
    tensors = ("sigma", "S", "kappae", "L", "PF", "RH")
    components = args.components
    if args.quantity in tensors:
        if args.components is None:
            lexit("{} is a tensor but no components have been specified"
                  .format(args.quantity))
    if args.quantity not in tensors:
        if components:
            warning(("{} is not a tensor. "
                     "The --components options will have no effect"
                     ).format(args.quantity))
        components = [()]
    elif args.quantity == "RH":
        if not all(i is None or len(i) == 3 for i in args.components):
            lexit("component specifications for the Hall tensor"
                  " need three indices")
    else:
        if not all(i is None or len(i) == 2 for i in args.components):
            lexit(("component specifications for {}"
                   " need two indices").format(args.quantity))
    # Prepare the abscissas, the third variable and the axis labels
    if args.abscissa == "T":
        x = Tr
        xlabel = r"$T\;\left[\mathrm{K}\right]$"
        z = mur[::args.subsample]
        zlabel0 = r"$\mu - \varepsilon_F = {:5g}\;\mathrm{{Ry}}$"
    else:
        x = mur - fermi
        xlabel = r"$\mu - \varepsilon_F\;\left[\mathrm{Ha}\right]$"
        z = Tr[::args.subsample]
        zlabel0 = r"$T = {:5g}\;\mathrm{{K}}$"
    ylabel = dict(
        cv=r"$c_v\;\left[\mathrm{J\,mol^{-1}\,K^{-1}}\right]$",
        n=r"$n\;\left[\mathrm{\left|e\right|\,uc^{-1}}\right]$",
        DOS=r"$\mathrm{DOS}\;\left[\mathrm{uc^{-1}}\right]$",
        sigma=r"$\sigma^{{\left({}\right)}}/\tau_0\;"
        r"\left[\mathrm{{ohm^{{-1}}\,m^{{-1}}\,"
        r"s^{{-1}}}}\right]$",
        S=r"$S^{{\left({}\right)}}\;" +
        r"\left[\mathrm{{V\,K^{{-1}}}}\right]$",
        kappae=r"$\kappa_e^{{\left({}\right)}}/\tau_0\;"
        r"\left[\mathrm{{W\,m^{{-1}}\,K^{{-1}}\,s^{{-1}}}}\right]$",
        L=r"$L^{{\left({}\right)}}\;"
        r"\left[\mathrm{{W\,\Omega\,K^{{-2}}}}\right]$",
        PF=r"$\left(S^2\sigma\right)^{{\left({}\right)}}/ \tau_0"
        r"\;\left[\mathrm{{W\,m^{{-1}}\,K^{{-2}}\,s^{{-1}}}}\right]$",
        RH=r"$R_H^{{\left({}\right)}}\;" +
        r"\left[\mathrm{{m^3\,C^{{-1}}}}\right]$")
    ylabel0 = ylabel[args.quantity]
    # Prepare the ordinate
    if args.quantity == "cv":
        y = cv * AVOGADRO / nformulas
    elif args.quantity == "n":
        y = N + data.nelect
    elif args.quantity == "DOS":
        y = sdos
    elif args.quantity == "sigma":
        y = cond
    elif args.quantity == "S":
        y = seebeck
    elif args.quantity == "kappae":
        y = kappa
    elif args.quantity == "L":
        L0 = 2.44e-8
        y = np.empty_like(cond)
        for iT in range(y.shape[0]):
            for imu in range(y.shape[1]):
                y[iT, imu] = la.solve(cond[iT, imu].T,
                                      kappa[iT, imu].T).T / Tr[iT]
    elif args.quantity == "PF":
        y = np.empty_like(cond)
        for iT in range(y.shape[0]):
            for imu in range(y.shape[1]):
                y[iT, imu] = (cond[iT, imu] @ seebeck[iT, imu]
                              @ seebeck[iT, imu])
    elif args.quantity == "RH":
        y = hall
    if args.abscissa == "T":
        y = y[:, ::args.subsample]
    else:
        y = y[::args.subsample, :]
    # Create the plots
    for c in components:
        if args.quantity in tensors:
            desc = "scalar" if c is None else "".join("xyz" [i] for i in c)
            ylabel = ylabel0.format(desc)
        else:
            ylabel = ylabel0
        plt.figure()
        for iz, zv in enumerate(z):
            zlabel = zlabel0.format(zv)
            if args.abscissa == "T":
                indices = [slice(None, None, None), iz]
            else:
                indices = [iz, slice(None, None, None)]
            if c is None:
                thisy = np.zeros_like(x)
                if args.quantity == "RH":
                    complements = ((0, 1, 2), (1, 2, 0), (2, 0, 1))
                else:
                    complements = [[i, i] for i in range(3)]
                for i in complements:
                    tindices = tuple(indices + list(i))
                    thisy += y[tindices]
                thisy /= float(len(complements))
            else:
                indices += c
                indices = tuple(indices)
                thisy = y[indices]
            plt.plot(x, thisy, label=zlabel)
        if args.abscissa == "mu":
            plt.axvline(x=0., color=PSEUDO_BLACK, lw=2)
        if args.quantity == "L":
            plt.axhline(y=L0, color=PSEUDO_BLACK, lw=2)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend(loc="best").draggable(True)
        plt.tight_layout()
    plt.show()


def parse_describe(args):
    """Parse command-line arguments for the describe subcommand."""
    data, metadata = BoltzTraP2.serialization.load_basic(args.file)
    print("Compound:", data.atoms.get_chemical_formula())
    print("Metadata:")
    for k in metadata:
        print("\t - {}: {}".format(k, metadata[k]))


def parse_fermisurface(args):
    """Parse command-line arguments for the fermisurface subcommand."""
    global terminal_colors
    # Try and load the data from the interpolation step
    data, equivalences, coeffs, metadata = (
        BoltzTraP2.serialization.load_calculation(args.bt2_file))
    info("sucessfully loaded " + args.bt2_file)
    lattvec = data.get_lattvec()
    # Rebuild the bands
    with TimerContext() as timer:
        ebands = BoltzTraP2.fite.getBTPbands(equivalences, coeffs, lattvec,
                                             True, args.nworkers)[0]
        deltat = timer.get_deltat()
        info("rebuilding the bands took {:.3g} s".format(deltat))
    # Print out some useful information
    if terminal_colors:
        colors = dict(
            green=colorama.Fore.GREEN,
            red=colorama.Fore.RED,
            bold=colorama.Style.BRIGHT,
            normal=colorama.Style.RESET_ALL)
    else:
        colors = dict(green="", red="", bold="", normal="")
    print(colors["bold"] + "─" * 72 + colors["normal"])
    print("""{bold}Keyboard and mouse interface:{normal}
{green}r{normal}: reset camera
{green}w{normal}: switch between surface and wireframe mode
{green}e{normal}: exit

{red}left mouse button{normal}: rotate around in-screen axes
{green}CTRL{normal} + {red}left mouse button{normal}: rotate around through-screen axis
{red}right mouse button{normal}: zoom
{red}middle mouse button{normal}: move / pick a point""".format(**colors))
    print(colors["bold"] + "─" * 72 + colors["normal"])
    # Call the function in charge of setting up the representation.
    # The function will only return after the user closes the interface.
    BoltzTraP2.fermisurface.plot_fermisurface(data, equivalences, ebands,
                                              args.mu, args.nworkers)


def parse_arguments():
    """Parse the command-line arguments and call the right subcommand."""
    # Create the main parser and perform the general setup.
    parser = argparse.ArgumentParser(
        description="BoltzTraP2, a program for interpolating band structures "
        "and calculating semi-classical transport coefficients",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="increase the verbosity level")
    parser.add_argument(
        "-n",
        "--nworkers",
        type=positive_integer,
        default=1,
        help="number of processes to span for parallel operations")
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version="BoltzTraP2 {}".format(PROGRAM_VERSION),
        help="print version information")
    # Create the subcommand parsers.
    subparsers = parser.add_subparsers(
        title="available subcommands", metavar="command [options] ...")
    subparsers.required = True
    interpolate_parser = subparsers.add_parser(
        "interpolate",
        help="interpolate the DFT input and save the results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    integrate_parser = subparsers.add_parser(
        "integrate",
        help="compute transport coefficients based on the interpolated bands",
        epilog="The reference for the chemical potential is electroneutrality.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    plotbands_parser = subparsers.add_parser(
        "plotbands",
        help="create a plot of the interpolated band structure",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    plot_parser = subparsers.add_parser(
        "plot",
        help="create a plot of the integration results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    describe_parser = subparsers.add_parser(
        "describe",
        help="provide basic information about a .bt2 or a .btj file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The 'fermisurface' subcommand is only acailable if vtk can be imported
    if BoltzTraP2.fermisurface.available:
        fermisurface_parser = subparsers.add_parser(
            "fermisurface",
            help="create an interactive visualization of the Fermi surface",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Options for the "interpolate" subcommand.
    interpolate_parser.add_argument(
        "-o",
        "--output",
        default="interpolation.bt2",
        help="name of the output file")
    interpolate_parser.add_argument(
        "-d",
        "--derivatives",
        action="store_true",
        help="use the derivatives of the bands in the interpolation")
    interpolate_parser.add_argument(
        "-e",
        "--emin",
        type=float,
        default=-0.2,
        help="bands with points below this energy"
        " (relative to the Fermi level) will be ignored")
    interpolate_parser.add_argument(
        "-E",
        "--emax",
        type=float,
        default=+0.2,
        help="bands with points above this energy"
        " (relative to the Fermi level) will be ignored")
    interpolate_parser.add_argument(
        "-a",
        "--absolute",
        action="store_true",
        help="interpret -e/--emin and -E/--emax  as absolute energies")
    kpoints_group = interpolate_parser.add_mutually_exclusive_group(
        required=True)
    kpoints_group.add_argument(
        "-k",
        "--kpoints",
        type=positive_integer,
        help="approximate number of irreducible k points to be used in the"
        " interpolation")
    kpoints_group.add_argument(
        "-m",
        "--multiplier",
        type=positive_integer,
        help="enhancement factor for the number of irreducible k points in the"
        " interpolation vs the DFT input")
    interpolate_parser.add_argument("directory", help="path to the DFT inputs")
    interpolate_parser.set_defaults(func=parse_interpolate)
    # Options for the "integrate" subcommand.
    integrate_parser.add_argument(
        "bt2_file", help="saved output from an interpolation step")
    integrate_parser.add_argument(
        "temperature",
        type=temperature_argument,
        help="single temperature, 3-tuple (Tmin, Tmax, nT), or list of"
        " temperatures to explore [K]")
    integrate_parser.add_argument(
        "-p", "--prefix", help="prefix for the output files")
    integrate_parser.add_argument(
        "-b",
        "--bins",
        type=positive_integer,
        help="number of bins for the DOS")
    scattering_group = integrate_parser.add_mutually_exclusive_group()
    scattering_group.add_argument(
        "-t",
        "--tau",
        action="store_const",
        dest="scattering_model",
        const="uniform_tau",
        help="work under a uniform-relaxation-time approximation")
    scattering_group.add_argument(
        "-l",
        "--lambda",
        action="store_const",
        dest="scattering_model",
        const="uniform_lambda",
        help="work under a uniform-mean-free-path approximation")
    integrate_parser.set_defaults(scattering_model="tau", func=parse_integrate)
    # Options for the "plotbands" subcommand.
    plotbands_parser.add_argument(
        "bt2_file", help="saved output from an interpolation step")
    plotbands_parser.add_argument(
        "-k",
        "--nkpoints",
        type=positive_integer,
        default=50,
        help="number of k points to sample along each the path")
    plotbands_parser.add_argument(
        "kpath",
        help="path in reciprocal space [Python list sytax, "
        "segments separated by None]")
    plotbands_parser.set_defaults(func=parse_plotbands)
    # Options for the "plot" subcommand
    plot_parser.add_argument(
        "btj_file", help="saved output from an integration step")
    plot_parser.add_argument(
        "quantity",
        choices=("cv", "n", "DOS", "sigma", "S", "kappae", "L", "PF", "RH"),
        help="quantity to be plotted")
    plot_parser.add_argument(
        "-c",
        "--components",
        type=components_argument,
        help="components to be plotted in the case of tensor quantitites")
    abscissa_group = plot_parser.add_mutually_exclusive_group()
    abscissa_group.add_argument(
        "-u",
        "--mu",
        action="store_const",
        dest="abscissa",
        const="mu",
        help="use the chemical potential as the abscissa in the plots")
    abscissa_group.add_argument(
        "-T",
        "--temperature",
        action="store_const",
        dest="abscissa",
        const="T",
        help="use the temperature as the abscissa in the plots")
    plot_parser.add_argument(
        "-s",
        "--subsample",
        type=positive_integer,
        default=1,
        help="plot a curve only for one out of each SUBSAMPLE values of the"
        " variable not chosen as the abscissa")
    plot_parser.set_defaults(abscissa="mu", func=parse_plot)
    # Options for the "describe" subcommand
    describe_parser.add_argument("file", help=".bt2 or .btj file")
    describe_parser.set_defaults(func=parse_describe)
    # Options for the "fermisurface" subcommand, if present
    if BoltzTraP2.fermisurface.available:
        fermisurface_parser.add_argument(
            "bt2_file", help="saved output from an interpolation step")
        fermisurface_parser.add_argument(
            "mu", type=float, help="chemical potential [Ha]")
        fermisurface_parser.set_defaults(func=parse_fermisurface)
    # Parse the only global argument and delegate on the chosen parser.
    for i, arg in enumerate(sys.argv):
        # This is a quick workaround for Python's inconvenient policy of
        # parsing negative numbers as short arguments.
        if arg.startswith("-") and len(arg) > 1 and (arg[1].isdigit() or
                                                     arg[1] == "."):
            sys.argv[i] = " " + arg
    args = parser.parse_args()
    set_logging_level(args.verbose)
    args.func(args)


def btp2_main():
    """Entry point for the btp2 script."""
    # Set up logging for the program
    logging.basicConfig(
        level=logging.ERROR, format="{levelname:8s}│ {message:s}", style="{")
    # Each subcommand is handled by a specialized function and there is little
    # left to do here except calling the command-line parser.
    parse_arguments()
