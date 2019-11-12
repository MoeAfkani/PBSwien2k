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

import multiprocessing as mp
import multiprocessing.sharedctypes

import numpy as np
import scipy as sp
# Try to use the FFT implementation in PyFFTW if available, or fall back on
# numpy's.
try:
    import pyfftw.interfaces.numpy_fft as npf
except ModuleNotFoundError:
    from BoltzTraP2.misc import warning
    warning("you can install pyfftw to get better FFT performance")
    import numpy.fft as npf

from BoltzTraP2 import sphere
from BoltzTraP2.units import *


def exp_worker(equivalences, kp, lattvec, mommat, sphase, sphaseR, index):
    """Worker function implementing the inner loop of fitde3D.

    Args:
        equivalences: list of k-point equivalence classes in direct coordinates
        kp: irreducible set of k points
        lattvec: lattice vectors of the system
        mommat: momentum matrix elements if available, or None otherwise
        sphase: shared-memory array to fill with results
        sphaseR: shared-memory array to fill with results
        index: index into sphase and sphaseR

    Returns:
        None. The results of the calculation are put in sphase and sphaseR.
    """
    if mommat is not None:
        dimR = 4 * len(kp) - 1
    else:
        dimR = len(kp) - 1
    allphase = np.frombuffer(sphase, dtype=np.complex)
    allphase.shape = (len(kp), -1)
    allphaseR = np.frombuffer(sphaseR, dtype=np.complex)
    allphaseR.shape = (dimR, -1)

    tpii = 2j * np.pi
    nstar = np.array([len(equiv) for equiv in equivalences])

    phase = np.zeros((len(kp), len(equivalences)), dtype=complex)
    if mommat is not None:
        phaseR = np.zeros((3, len(kp), len(equivalences)), dtype=complex)
    for j, equiv in enumerate(equivalences):
        phase0 = np.exp(tpii * kp @ equiv.T)
        phase[:, j] = np.sum(phase0, axis=1)
        if mommat is not None:
            vv = 1j * (lattvec @ equiv.T)
            phaseR[:, :, j] = vv @ phase0.T
    phase /= nstar

    if mommat is not None:
        phaseR /= nstar
        phaseR = np.vstack((phase[:-1] - phase[-1], phaseR[0], phaseR[1],
                            phaseR[2]))
    else:
        phaseR = phase[:-1] - phase[-1]

    allphase[:, index:index + len(equivalences)] = phase
    allphaseR[:, index:index + len(equivalences)] = phaseR


def fitde3D(data, equivalences, nworkers=1):
    """Obtain the interpolation coefficients from DFT data.

    Args:
        data: DFTdata object containing all relevant input
        equivalences: list of k-point equivalence classes in direct coordinates
        nworkers: number of working processes to span

    Returns:
        A set of interpolation coefficients as an array whose first dimension
        runs over bands.
    """
    kp = data.kpoints
    ene = data.ebands
    mommat = data.mommat
    lattvec = data.get_lattvec()

    C1 = .75
    C2 = .75
    Rvec = np.array([equiv[0] for equiv in equivalences])

    R = np.linalg.norm(Rvec @ lattvec.T, axis=1)
    De = ene.T[:-1] - ene.T[-1]
    if mommat is not None:
        for i in range(3):
            De = np.vstack((De, mommat[:, :, i]))
    X2 = (R / R[1])**2
    rhoi = 1. / ((1. - C1 * X2)**2 + C2 * X2**3)
    rhoi[0] = 0.

    # Split the equivalence classes among worker processes
    chunk_size = int(np.ceil(len(equivalences) / float(nworkers)))
    workers = []
    # Use chunks of shared memory to store the results
    # Note that double-sized real arrays are used to transparently store
    # complex data
    if mommat is not None:
        dimR = 4 * len(kp) - 1
    else:
        dimR = len(kp) - 1
    sphase = mp.sharedctypes.RawArray("d", len(kp) * len(equivalences) * 2)
    phase = np.frombuffer(sphase, dtype=np.complex)
    phase.shape = (len(kp), len(equivalences))
    sphaseR = mp.sharedctypes.RawArray("d", dimR * len(equivalences) * 2)
    phaseR = np.frombuffer(sphaseR, dtype=np.complex)
    phaseR.shape = (dimR, len(equivalences))
    for i in range(0, len(equivalences), chunk_size):
        subequivalences = equivalences[i:i + chunk_size]
        workers.append(
            mp.Process(
                target=exp_worker,
                args=(subequivalences, kp, lattvec, mommat, phase, phaseR, i)))
    for w in workers:
        w.start()
    # Wait for the worker processes to finish their job
    for w in workers:
        w.join()
    # Put the results in arrays with the right layout before the final
    # operations.
    pRR = np.ascontiguousarray(phaseR[:, 1:].real)
    pRI = np.ascontiguousarray(phaseR[:, 1:].imag)
    qRR = (pRR * rhoi[1:]).T
    qRI = np.ascontiguousarray(qRR.imag)
    qRR = np.ascontiguousarray(qRR.real)
    Hmat = pRR @ qRR + pRI @ qRI
    rlambda = np.linalg.lstsq(Hmat, De)[0]

    coeffs = rhoi * (rlambda.T @ phaseR)
    coeffs[:, 0] = ene.T[-1] - coeffs[:, 1:] @ phase[-1, 1:]

    return coeffs


def fft_worker(equivalences, sallvec, dims, iqueue, oqueue, curvature=False):
    """Thin wrapper around FFTev and FFTc to be used as a worker function.

    Args:
        equivalences: list of k-point equivalence classes in direct coordinates
        sallvec: Cartesian coordinates of all k points as a 1D vector stored
                    in shared memory.
        dims: upper bound on the dimensions of the k-point grid
        iqueue: input multiprocessing.Queue used to read bad indices
            and coefficients.
        oqueue: output multiprocessing.Queue where all results of the
            interpolation are put. Each element of the queue is a 4-tuple
            of the form (index, eband, vvband, cband), containing the band
            index, the energies, the v x v outer product and the curvatures
            if requested.
        curvature: should the band curvature be calculated?

    Returns:
        None. The results of the calculation are put in oqueue.
    """
    iu0 = np.triu_indices(3)
    il1 = np.tril_indices(3, -1)
    iu1 = np.triu_indices(3, 1)
    allvec = np.frombuffer(sallvec)
    allvec.shape = (-1, 3)
    if curvature:
        invmass = np.zeros((3, 3, np.prod(dims)))
    while True:
        task = iqueue.get()
        if task is None:
            break
        else:
            index, bandcoeff = task
        eband, vb = FFTev(equivalences, bandcoeff, allvec, dims)
        vvband = np.zeros((3, 3, np.prod(dims)))
        vvband[iu0[0], iu0[1]] = vb[iu0[0]] * vb[iu0[1]]
        vvband[il1[0], il1[1]] = vvband[iu1[0], iu1[1]]
        if curvature:
            invmass[iu0] = FFTc(equivalences, bandcoeff, allvec, dims)
            invmass[il1] = invmass[iu1]
            cband = np.zeros((3, 3, 3, np.prod(dims)))
            for i in range(3):
                for j in range(3):
                    cband[i, j] = np.cross(vvband[i].T, invmass[j].T).T
        else:
            cband = None
        oqueue.put((index, eband, vvband, cband))


def getBTPbands(equivalences, coeffs, lattvec, curvature=False, nworkers=1):
    """Rebuild the full energy bands from the interpolation coefficients.

    Args:
        equivalences: list of k-point equivalence classes in direct coordinates
        coeffs: interpolation coefficients
        lattvec: lattice vectors of the system
        curvature: should the band curvature be calculated?
        nworkers: number of working processes to span

    Returns:
        A 3-tuple (eband, vvband, cband): energy bands, v x v outer product
        of the velocities, and curvature of the bands (if requested). The
        shapes of those arrays are (nbands, nkpoints), (nbands, 3, 3, nkpoints)
        and (nbands, 3, 3, 3, nkpoints), where nkpoints is the total number of
        k points on the grid. If curvature is None, so will the third element
        of the tuple.
    """
    dallvec = np.vstack(equivalences)
    sallvec = mp.sharedctypes.RawArray("d", dallvec.shape[0] * 3)
    allvec = np.frombuffer(sallvec)
    allvec.shape = (-1, 3)
    dims = 2 * np.max(np.abs(dallvec), axis=0) + 1
    np.matmul(dallvec, lattvec.T, out=allvec)
    eband = np.zeros((len(coeffs), np.prod(dims)))
    vvband = np.zeros((len(coeffs), 3, 3, np.prod(dims)))
    if curvature:
        cband = np.zeros((len(coeffs), 3, 3, 3, np.prod(dims)))
    else:
        cband = None

    # Span as many worker processes as needed, put all the bands in the queue,
    # and let them work until all the required FFTs have been computed.
    workers = []
    iqueue = mp.Queue()
    oqueue = mp.Queue()
    for iband, bandcoeff in enumerate(coeffs):
        iqueue.put((iband, bandcoeff))
    # The "None"s at the end of the queue signal the workers that there are
    # no more jobs left and they must therefore exit.
    for i in range(nworkers):
        iqueue.put(None)
    for i in range(nworkers):
        workers.append(
            mp.Process(
                target=fft_worker,
                args=(equivalences, sallvec, dims, iqueue, oqueue, curvature)))
    for w in workers:
        w.start()
    # The results of the FFTs are processed as soon as they are ready.
    for r in range(len(coeffs)):
        iband, eband[iband], vvband[iband], cb = oqueue.get()
        if curvature:
            cband[iband] = cb
    for w in workers:
        w.join()
    if cband is not None:
        cband = cband.real
    return eband.real, vvband.real, cband


def FFTev(equivalences, bandcoeff, allvec, dims):
    """Rebuild a single band from the interpolation coefficients.

    Args:
        equivalences: list of k-point equivalence classes in direct coordinates
        bandcoeff: interpolation coefficients for the band
        allvec: Cartesian coordinates of all k points
        dims: upper bound on the dimensions of the k-point grid

    Returns:
        A 2-tuple (eb, vg) with the energies and group velocities of the
        reconstructed bands.
    """
    egrid = np.zeros(dims, dtype=complex)
    vgrid = np.zeros([3] + list(dims), dtype=complex)
    i = 0
    for coeff, equiv in zip(bandcoeff, equivalences):
        j = i + len(equiv)
        c = coeff / len(equiv)
        egrid[equiv[:, 0], equiv[:, 1], equiv[:, 2]] = c
        vgrid[:, equiv[:, 0], equiv[:, 1], equiv[:, 2]] = allvec[i:j, :].T * c
        i = j
    vgrid *= 1j
    eb = np.prod(dims) * npf.ifftn(egrid).real.flatten()
    vg = [
        np.prod(dims) * npf.ifftn(vgrid[0]).real.flatten(),
        np.prod(dims) * npf.ifftn(vgrid[1]).real.flatten(),
        np.prod(dims) * npf.ifftn(vgrid[2]).real.flatten()
    ]
    return eb, np.array(vg)


def FFTc(equivalences, bandcoeff, allvec, dims):
    """Rebuild the curvature of a band from the interpolation coefficients.

    Args:
        equivalences: list of k-point equivalence classes in direct coordinates
        bandcoeff: interpolation coefficients for the band
        allvec: Cartesian coordinates of all k points
        dims: upper bound on the dimensions of the k-point grid

    Returns:
        An array used to rebuild the curvatures in fft_worker.
    """
    cgrid = np.zeros([6] + list(dims), dtype=complex)
    i = 0
    ii1, ii2 = np.triu_indices(3)
    for coeff, equiv in zip(bandcoeff, equivalences):
        j = i + len(equiv)
        c = coeff / len(equiv)
        eq = allvec[i:j, :]
        cgrid[:, equiv[:, 0], equiv[:, 1], equiv[:, 2]] = -(eq[:, ii1] *
                                                            eq[:, ii2]).T * c
        i = j
    cg = [npf.ifftn(cgrid[i]).real.flatten() for i in range(6)]
    return np.prod(dims) * np.array(cg)


def getBands(kp, equivalences, lattvec, coeffs):
    """Sample the energy bands at particular points from the interpolation
    coefficients.

    Args:
        kp: list of k points on which to sample the bands
        equivalences: list of k-point equivalence classes in direct coordinates
            used for the interpolation.
        lattvec: lattice vectors of the system
        coeffs: interpolation coefficients

    Returns:
        A 2-tuple (eband, vband): energy bands and group velocities
        reconstructed from the interpolation. Theshapes of those arrays are
        (nbands, nkpoints), (nbands, nkpoints) and (3, nbands, nkpoints).
    """
    tpii = 2j * np.pi
    phase = np.zeros((len(equivalences), len(kp)), dtype=complex)
    phaseR = np.zeros((len(equivalences), len(kp), 3), dtype=complex)
    for j, equiv in enumerate(equivalences):
        phase0 = np.exp(tpii * kp @ equiv.T)
        phase[j] = np.sum(phase0, axis=1)
        vv = 1j * (lattvec @ equiv.T)
        phaseR[j] = phase0 @ vv.T
    nstar = np.array([len(equiv) for equiv in equivalences])
    phase = (phase.T / nstar).T
    phaseR = (phaseR.T / nstar).T
    egrid = coeffs @ phase
    vgrid = []
    for i in range(3):
        vgrid += [coeffs @ phaseR[:, :, i]]
    return egrid.real, np.array(vgrid).real
