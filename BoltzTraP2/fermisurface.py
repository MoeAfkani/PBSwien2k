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

import itertools
import logging

import numpy as np
import scipy as sp
import scipy.linalg as la
import scipy.spatial
import matplotlib.colors as colors

from BoltzTraP2.misc import info
from BoltzTraP2.misc import warning
from BoltzTraP2.misc import TimerContext

# If VTK cannot be impoorted, the module-level variable "available" will be set
# to False and the features of this submodule will fail to work.
try:
    import vtk
    available = True
except ModuleNotFoundError:
    warning("vtk not found. The 'fermisurface' command will be disabled")
    available = False


def plot_fermisurface(data, equivalences, ebands, mu, nworkers=1):
    """Launch an interactive VTK representation of the Fermi surface.

    Make sure to check the module-level variable "available" before calling
    this function.

    Args:
        data: a DFTData object
        equivalences: list of k-point equivalence classes
        ebands: eband: (nbands, nkpoints) array with the band energies
        mu: initial value of the energy at which the surface will be plotted,
            with respect to data.fermi. This can later be changed by the user
            using an interactive slider.
        nworkers: number of worker processes to use for the band reconstruction

    Returns:
        None.
    """
    lattvec = data.get_lattvec()
    rlattvec = 2. * np.pi * la.inv(lattvec).T

    # Obtain the first Brillouin zone as the Voronoi polyhedron of Gamma
    points = []
    for ijk0 in itertools.product(range(5), repeat=3):
        ijk = [i if i <= 2 else i - 5 for i in ijk0]
        points.append(rlattvec @ np.array(ijk))
    voronoi = scipy.spatial.Voronoi(points)
    region_index = voronoi.point_region[0]
    vertex_indices = voronoi.regions[region_index]
    vertices = voronoi.vertices[vertex_indices, :]
    # Compute a center and an outward-pointing normal for each of the facets
    # of the BZ
    facets = []
    for ridge in voronoi.ridge_vertices:
        if all(i in vertex_indices for i in ridge):
            facets.append(ridge)
    centers = []
    normals = []
    for f in facets:
        corners = np.array([voronoi.vertices[i, :] for i in f])
        center = corners.mean(axis=0)
        v1 = corners[0, :]
        for i in range(1, corners.shape[0]):
            v2 = corners[i, :]
            prod = np.cross(v1 - center, v2 - center)
            if not np.allclose(prod, 0.):
                break
        if np.dot(center, prod) < 0.:
            prod = -prod
        centers.append(center)
        normals.append(prod)

    # Get the extent of the regular grid in reciprocal space
    hdims = np.max(np.abs(np.vstack(equivalences)), axis=0)
    dims = 2 * hdims + 1

    class PointPickerInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):
        """Custom interaction style enabling the user to pick points on
        the screen.
        """

        def __init__(self, parent=None):
            """Simple constructor that adds an observer to the middle mouse
            button press.
            """
            self.AddObserver("MiddleButtonPressEvent", self.pick_point)

        def pick_point(self, obj, event):
            """Get the coordinates of the point selected with the middle mouse
            button, find the nearest data point, and print its direct
            coordinates.
            """
            interactor = self.GetInteractor()
            picker = interactor.GetPicker()
            pos = interactor.GetEventPosition()
            picker.Pick(
                pos[0], pos[1], 0,
                interactor.GetRenderWindow().GetRenderers().GetFirstRenderer())
            picked = np.array(picker.GetPickPosition())
            # Move the sphere to the new coordinates and make it visible
            sphere.SetCenter(*picked.tolist())
            sphere_actor.VisibilityOn()
            picked = la.solve(rlattvec, picked)
            print("Point picked:", picked)
            self.OnMiddleButtonDown()

    # Create the VTK representation of the grid
    sgrid = vtk.vtkStructuredGrid()
    sgrid.SetDimensions(*dims)
    spoints = vtk.vtkPoints()
    for ijk0 in itertools.product(*(range(0, d) for d in dims)):
        ijk = [
            ijk0[i] if ijk0[i] <= hdims[i] else ijk0[i] - dims[i]
            for i in range(len(dims))
        ]
        abc = np.array(ijk, dtype=np.float64) / np.array(dims)
        xyz = rlattvec @ abc
        spoints.InsertNextPoint(*xyz.tolist())
    sgrid.SetPoints(spoints)

    # Find the shortest distance between points to compute a good
    # radius for the selector sphere later.
    dmin = np.infty
    for i in range(3):
        abc = np.zeros(3)
        abc[i] = 1. / dims[i]
        xyz = rlattvec @ abc
        dmin = min(dmin, la.norm(xyz))

    ebands -= data.fermi
    emax = ebands.max(axis=1)
    emin = ebands.min(axis=1)
    # Remove all points outside the BZ
    for i, ijk0 in enumerate(itertools.product(*(range(0, d) for d in dims))):
        ijk = [
            ijk0[j] if ijk0[j] <= hdims[j] else ijk0[j] - dims[j]
            for j in range(len(dims))
        ]
        abc = np.array(ijk, dtype=np.float64) / np.array(dims)
        xyz = rlattvec @ abc
        for c, n in zip(centers, normals):
            if np.dot(xyz - c, n) > 0.:
                ebands[:, i] = np.nan
                break

    # Create a 2D chemical potential slider
    slider = vtk.vtkSliderRepresentation2D()
    slider.SetMinimumValue(emin.min())
    slider.SetMaximumValue(emax.max())
    slider.SetValue(mu)
    slider.SetTitleText("Chemical potential")
    slider.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider.GetPoint1Coordinate().SetValue(0.1, 0.9)
    slider.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider.GetPoint2Coordinate().SetValue(0.9, 0.9)
    slider.GetTubeProperty().SetColor(*colors.hex2color("#2e3436"))
    slider.GetSliderProperty().SetColor(*colors.hex2color("#a40000"))
    slider.GetCapProperty().SetColor(*colors.hex2color("#babdb6"))
    slider.GetSelectedProperty().SetColor(*colors.hex2color("#a40000"))
    slider.GetTitleProperty().SetColor(*colors.hex2color("#2e3436"))

    # Find all the isosurfaces with energy equal to the threshold
    allcontours = []
    with TimerContext() as timer:
        fermiactors = []
        for band in ebands:
            sgridp = vtk.vtkStructuredGrid()
            sgridp.DeepCopy(sgrid)
            # Feed the energies to VTK
            scalar = vtk.vtkFloatArray()
            for i in band:
                scalar.InsertNextValue(i)
            sgridp.GetPointData().SetScalars(scalar)
            # Estimate the isosurfaces
            contours = vtk.vtkMarchingContourFilter()
            contours.SetInputData(sgridp)
            contours.UseScalarTreeOn()
            contours.SetValue(0, mu)
            contours.ComputeNormalsOff()
            contours.ComputeGradientsOff()
            allcontours.append(contours)

            # Use vtkStrippers to plot the surfaces faster
            stripper = vtk.vtkStripper()
            stripper.SetInputConnection(contours.GetOutputPort())

            # Compute the normals to the surfaces to obtain better lighting
            normals = vtk.vtkPolyDataNormals()
            normals.SetInputConnection(stripper.GetOutputPort())
            normals.ComputeCellNormalsOn()
            normals.ComputePointNormalsOn()

            # Create a mapper and an actor for the surfaces
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(normals.GetOutputPort())
            mapper.ScalarVisibilityOff()
            fermiactors.append(vtk.vtkActor())
            fermiactors[-1].SetMapper(mapper)
            fermiactors[-1].GetProperty().SetColor(*colors.hex2color(
                "#a40000"))
        deltat = timer.get_deltat()
        info("building the surfaces took {:.3g} s".format(deltat))

    # Represent the BZ as a polyhedron in VTK
    points = vtk.vtkPoints()
    for v in voronoi.vertices:
        points.InsertNextPoint(*v)
    fids = vtk.vtkIdList()
    fids.InsertNextId(len(facets))
    for f in facets:
        fids.InsertNextId(len(f))
        for i in f:
            fids.InsertNextId(i)
    fgrid = vtk.vtkUnstructuredGrid()
    fgrid.SetPoints(points)
    fgrid.InsertNextCell(vtk.VTK_POLYHEDRON, fids)

    # Create an actor and a mapper for the BZ
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(fgrid)
    bzactor = vtk.vtkActor()
    bzactor.SetMapper(mapper)
    bzactor.GetProperty().SetColor(*colors.hex2color("#204a87"))
    bzactor.GetProperty().SetOpacity(0.2)

    # Create a visual representation of the selected point, and hide
    # it for the time being.
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(dmin / 2.)
    sphere_mapper = vtk.vtkPolyDataMapper()
    sphere_mapper.SetInputConnection(sphere.GetOutputPort())
    sphere_mapper.ScalarVisibilityOff()
    sphere_actor = vtk.vtkActor()
    sphere_actor.SetMapper(sphere_mapper)
    sphere_actor.GetProperty().SetColor(*colors.hex2color("#f57900"))
    sphere_actor.VisibilityOff()

    # Create a VTK window and other elements of an interactive scene
    renderer = vtk.vtkRenderer()
    renderer.AddActor(bzactor)
    renderer.AddActor(sphere_actor)
    for f in fermiactors:
        renderer.AddActor(f)
    renderer.ResetCamera()
    renderer.GetActiveCamera().Zoom(5.)
    renderer.SetBackground(1., 1., 1.)

    window = vtk.vtkRenderWindow()
    window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetInteractorStyle(PointPickerInteractorStyle())
    interactor.SetRenderWindow(window)

    # Add a set of axes
    axes = vtk.vtkAxesActor()
    assembly = vtk.vtkPropAssembly()
    assembly.AddPart(axes)
    marker = vtk.vtkOrientationMarkerWidget()
    marker.SetOrientationMarker(assembly)
    marker.SetInteractor(interactor)
    marker.SetEnabled(1)
    marker.InteractiveOff()

    def callback(obj, ev):
        """Update the isosurface with a new value"""
        mu = obj.GetRepresentation().GetValue()
        for e, E, c, a in zip(emin, emax, allcontours, fermiactors):
            visible = e <= mu and E >= mu
            a.SetVisibility(visible)
            if visible:
                c.SetValue(0, mu)

    # Add the slider widget
    widget = vtk.vtkSliderWidget()
    widget.SetInteractor(interactor)
    widget.SetRepresentation(slider)
    widget.SetAnimationModeToJump()
    widget.EnabledOn()
    widget.AddObserver(vtk.vtkCommand.InteractionEvent, callback)

    # Launch the visualization
    interactor.Initialize()
    window.Render()
    interactor.Start()
