#!/usr/bin/env python

import time
import math
import itertools
import numpy as np
import scipy.io as sio
from scipy import signal
import vtk
from vtk.util.numpy_support import numpy_to_vtk

def makeStiffnessMatrixUniqueElements(material_properties):
    lam = material_properties["lame_1"]
    mu = material_properties["lame_2"]
    
    C = np.zeros((21))
    C[0] = lam+2*mu
    C[1] = lam
    C[2] = lam
    C[6] = lam+2*mu
    C[7] = lam
    C[11] = lam+2*mu
    C[15] = mu
    C[18] = mu
    C[20] = mu
    return C

def makeAngleBracket(bracket_geometry, material_properties, timing, water_column):

    ts1 = time.time()
    ts = ts1

    dimensions = [bracket_geometry["bracket_length"], bracket_geometry["bracket_height"], bracket_geometry["bracket_height"]]
    ds = bracket_geometry["cell_dimensions"]
    bounding_dimensions_cells = []
    for i in range(3):
        bounding_dimensions_cells.append(int(np.ceil(dimensions[i]/ds[i])))
        pass

    geometry = np.zeros((bounding_dimensions_cells[0]+2, bounding_dimensions_cells[1]+2, bounding_dimensions_cells[2]+2), dtype=np.uint8)
    thickness_cells_y = int(np.ceil(bracket_geometry["bracket_thickness"]/ds[1]))
    thickness_cells_z = int(np.ceil(bracket_geometry["bracket_thickness"]/ds[2]))
    geometry[1:-1,1:1+thickness_cells_z,1:-1] = 1
    geometry[1:-1,1:-1,1:1+thickness_cells_z] = 1

    print("Geometry volume is %d x %d x %d cells." % (bounding_dimensions_cells[0], bounding_dimensions_cells[1], bounding_dimensions_cells[2]))
    print("Geometry volume is %f x %f x %f mm with an error of %f, %f, and  %f. mm" % (bounding_dimensions_cells[0]*ds[0]*1000, bounding_dimensions_cells[1]*ds[1]*1000, bounding_dimensions_cells[2]*ds[2]*1000, (bounding_dimensions_cells[0]*ds[0] - dimensions[0])*1000, (bounding_dimensions_cells[1]*ds[1] - dimensions[1])*1000, (bounding_dimensions_cells[2]*ds[2] - dimensions[2])*1000))

    tf = time.time()
    print("Geometry Made in %f seconds!" % (tf-ts))

    ## Create VTK Geometry and Associate Material Properties
    stiffness_matrix_values = makeStiffnessMatrixUniqueElements(material_properties)
    rho = material_properties["density"]
    number_of_cells = np.count_nonzero(geometry)
    rho_arr = np.ones(number_of_cells)
    C_arr = np.ones((21, number_of_cells))
    points = vtk.vtkPoints()
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.Allocate(number_of_cells)
    pts = np.ones((bounding_dimensions_cells[0]+1, bounding_dimensions_cells[1]+1, bounding_dimensions_cells[2]+1))*-1
    pindx = 0
    cell_idx = 0
    for i in range(bounding_dimensions_cells[0]+1):
        for j, k in itertools.product(range(bounding_dimensions_cells[1]+1), range(bounding_dimensions_cells[2]+1)):
            test = np.sum(geometry[i:i+2, j:j+2, k:k+2])
            if test >= 1:
                pts[i, j, k] = pindx
                pindx = pindx + 1
                points.InsertPoint(int(pts[i, j, k]), [i*ds[0], j*ds[1], k*ds[2]])
                if geometry[i, j, k] > 0:
                    cube = [int(pts[i-1, j-1, k-1]), int(pts[i, j-1, k-1]), int(pts[i, j, k-1]), int(pts[i-1, j, k-1]), int(pts[i-1, j-1, k]), int(pts[i, j-1, k]), int(pts[i, j, k]), int(pts[i-1, j, k])]
                    # [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1],[0, 1, 1]]
                    if (np.array(cube) < 0).any():
                        print("cube failure at %i,%i,%i\n"%(i, j, k))
                        pass
                    ugrid.InsertNextCell(vtk.VTK_HEXAHEDRON, 8, cube)
                    C_arr[:, cell_idx] = stiffness_matrix_values
                    rho_arr[cell_idx] = rho
                    cell_idx = cell_idx + 1
                    pass
                pass
            pass
        print("Slice %d of %d done!" % (i, bounding_dimensions_cells[0]+1))
        pass
    ugrid.SetPoints(points)
    ts = tf
    tf = time.time()
    print("Unstructured Grid Made in %f seconds!" % (tf-ts))

    array = numpy_to_vtk(rho_arr, deep=True)
    name = "density"
    array.SetName(name)
    ugrid.GetCellData().AddArray(array)

    C_idx = 0
    for i, j in itertools.combinations_with_replacement(range(6), 2):
        name = "C%d%d" %(i+1, j+1)
        array = numpy_to_vtk(C_arr[C_idx, :], deep=True)
        array.SetName(name)
        ugrid.GetCellData().AddArray(array)
        C_idx = C_idx +1
        pass


    ## Associate Timing Information
    dt = timing["dt"]
    tVec = np.arange(0, timing["simulation_duration"], dt)
    NT = int(tVec.shape[0]/timing["write_interval"])
    wo_arr = np.arange(NT+1)*timing["write_interval"]*dt

    array = numpy_to_vtk(np.array([dt]), deep=True)
    name = "dt"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(wo_arr, deep=True)
    name = "write_times"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)


    ## Define and Associate Transducers
    pulseLen = (water_column["cycle_count"]/water_column["center_frequency"]) # Pulse length (seconds) given # of cycles
    pulseSteps = int(np.ceil(pulseLen/dt))
    phase_rads = water_column["phase"]*np.pi/180 

    wf_tVec = np.arange(0, (pulseSteps)*dt, dt)
    wf = (np.sin(wf_tVec*water_column["center_frequency"]*2*np.pi+phase_rads)).astype(np.double) 
    window = signal.windows.hann(pulseSteps)
    wf = np.multiply(wf, window) 
    wf = wf/wf.max()

    array = numpy_to_vtk(np.array([1]), deep=True)
    name = "NWaterColumns"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(np.array([water_column["location"][0]]), deep=True)
    name = "WC0/XCenter"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(np.array([water_column["location"][1]]), deep=True)
    name = "WC0/YCenter"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(np.array([water_column["location"][2]]), deep=True)
    name = "WC0/ZCenter"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(np.array([water_column["radius"]]), deep=True)
    name = "WC0/Radius"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(wf_tVec, deep=True)
    name = "WC0/SignalTimes"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(wf, deep=True)
    name = "WC0/SignalSigmaXX"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(wf, deep=True)
    name = "WC0/SignalSigmaYY"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(wf, deep=True)
    name = "WC0/SignalSigmaZZ"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    ts = tf
    tf = time.time()
    print("Data Fields made in %f seconds!" % (tf-ts))

    ## Writing to File
    writer = vtk.vtkXMLDataSetWriter()
    writer.SetFileName("AL_Angle_Bracket_Example.vtu")
    writer.SetInputData(ugrid)
    writer.Write()
    ts = tf
    tf = time.time()
    print("File Written in %f seconds!" % (tf-ts))

    return

def main():
    #dimentions in m, s, Pa, kg, degrees, except write_interval which is in steps
    bracket_geometry = {
        "bracket_length": 50.8e-3, 
        "bracket_height": 25.4e-3, 
        "bracket_thickness": 3.175e-3, 
        "cell_dimensions": np.array([9.921875e-05, 9.921875e-05, 9.921875e-05])}

    material_properties = {
        "density": 2790, 
        "lame_1": 114.5e9,
        "lame_2":42.5e9,}

    timing = {
        "dt": 9.92e-9, 
        "write_interval": 50, 
        "simulation_duration": 99.734e-6}

    water_column =  {
        "radius": 3.96875, 
        "center_frequency": 200e3, 
        "cycle_count": 3, 
        "phase": 0, 
        "location": np.array([12.675e-3, 12.7e-3, 0.0])}
    
    makeAngleBracket(bracket_geometry, material_properties, timing, water_column)
    return


if __name__ == '__main__':
    main()