#!/usr/bin/env python

import time
import math
import itertools
import numpy as np
import scipy.io as sio
from scipy import signal
import vtk
from vtk.util.numpy_support import numpy_to_vtk

def ijTable(i,j):
    ij_table = np.array([[0,5,4],
                         [5,1,3],
                         [4,3,2]])
    return ij_table[i,j]

def getMatrixFromTriangle(triangle_elements_list):
    matrix_size = 0
    list_length = len(triangle_elements_list)
    while list_length > 0:
        matrix_size = matrix_size + 1
        list_length = list_length - matrix_size
        pass
    if list_length < 0:
        print("Wrong number of elements to define the upper right triangle of a square matrix (ex. 3, 6, 10, 15, 21 etc.)")
        return np.array([-1])

    matrix = np.zeros((matrix_size, matrix_size))
    index = 0
    for i, j in itertools.combinations_with_replacement(range(matrix.shape[0]), 2):
        matrix[i,j] = triangle_elements_list[index]
        index = index + 1
        pass
    return matrix

def getTriangleFromMatrix(matrix):
    triangle_elements_list = []
    if(matrix.shape[0]!=matrix.shape[1]):
        print("Not a square matrix!")
        return [-1]
    if(len(matrix.shape)>2):
        print("matrix has more than 2 dimensions!")
        return [-1]
    for i, j in itertools.combinations_with_replacement(range(matrix.shape[0]), 2):
        triangle_elements_list.append(matrix[i][j])
        pass
    return triangle_elements_list

def ijklFromIj(matrix_ij):
    if(matrix_ij.shape != (6,6)):
        print("matrix must be 6x6")
        return -1
    matrix_ijkl = np.zeros((3,3,3,3))
    for i, j, k, l in itertools.product(range(3),repeat=4):
        matrix_ijkl[i,j,k,l] = matrix_ij[ijTable(i,j),ijTable(k,l)]
        pass
    return matrix_ijkl

def ijFromIjkl(matrix_ijkl):
    if(matrix_ijkl.shape != (3,3,3,3)):
        print("matrix must be 3x3x3x3")
        return -1
    matrix_ij = np.zeros((6,6))
    for i, j, k, l in itertools.product(range(3),repeat=4):
        matrix_ij[ijTable(i,j),ijTable(k,l)] = matrix_ijkl[i,j,k,l]
        pass
    return matrix_ij

def rotateStiffnessMatrix(original_stiffness_matrix, rotation):
    ijkl_old = ijklFromIj(original_stiffness_matrix)
    alpha = rotation*np.pi/180
    R = np.identity(3)
    R[0,0] = np.cos(alpha)
    R[0,1] = -1*np.sin(alpha)
    R[1,1] = np.cos(alpha)
    R[1,0] = np.sin(alpha)
    ijkl_new = np.zeros((3,3,3,3))
    for i, j, k, l, x, y, z, w in itertools.product(range(3),repeat=8):
        ijkl_new[i,j,k,l] = ijkl_new[i,j,k,l] + R[i,x]*R[j,y]*R[k,z]*R[l,w]*ijkl_old[x,y,z,w]
        pass
    new_stiffness_matrix = ijFromIjkl(ijkl_new)
    return new_stiffness_matrix

def makeLaminate(laminate_geometry, material_properties, timing, transducer):

    ts1 = time.time()
    ts = ts1

    layup = laminate_geometry["ply_orientations"]
    number_of_plies = len(layup)
    dimensions = [laminate_geometry["lamiant_dimensions"][0], laminate_geometry["lamiant_dimensions"][1], laminate_geometry["ply_thickness"]*number_of_plies]
    ds = laminate_geometry["cell_dimensions"]
    bounding_dimensions_cells = []
    for i in range(3):
        bounding_dimensions_cells.append(int(np.ceil(dimensions[i]/ds[i])))
        pass
    print("Geometry volume is %d x %d x %d cells." % (bounding_dimensions_cells[0], bounding_dimensions_cells[1], bounding_dimensions_cells[2]))
    print("Geometry volume is %f x %f x %f mm with an error of %f, %f, and  %f. mm" % (bounding_dimensions_cells[0]*ds[0]*1000, bounding_dimensions_cells[1]*ds[1]*1000, bounding_dimensions_cells[2]*ds[2]*1000, (bounding_dimensions_cells[0]*ds[0] - dimensions[0])*1000, (bounding_dimensions_cells[1]*ds[1] - dimensions[1])*1000, (bounding_dimensions_cells[2]*ds[2] - dimensions[2])*1000))

    geometry = np.zeros((bounding_dimensions_cells[0]+2, bounding_dimensions_cells[1]+2, bounding_dimensions_cells[2]+2), dtype=np.uint8)
    cells_per_ply = int(np.ceil(laminate_geometry["ply_thickness"]/ds[2]))
    C = material_properties["lamina_stiffness_matrix"]
    C_list = []
    rho_list = []
    for ply_i in range(number_of_plies):
        C_list.append(rotateStiffnessMatrix(C, layup[ply_i] ))
        rho_list.append(material_properties["density"])
        for i in range(cells_per_ply):
            material_id = ply_i+1
            geometry[1:-1, 1:-1, ply_i*cells_per_ply+i+1] = material_id
            pass
        print("Ply %d material properties calculated," % (ply_i+1))
        pass

    ### ADD a DELAMINATION ####
    # geometry[100:150,100:150,17]=0


    tf = time.time()
    print("Geometry Made in %f seconds!" % (tf-ts))

    ## Create VTK Geometry and Associate Material Properties
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
                    C_arr[:, cell_idx] = getTriangleFromMatrix(C_list[geometry[i, j, k]-1])
                    rho_arr[cell_idx] = rho_list[geometry[i, j, k]-1]
                    cell_idx = cell_idx + 1
                    pass
                pass
            pass
        print("Slice %d of %d done!" % (i+1, bounding_dimensions_cells[0]+1))
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
    pulseLen = (transducer["cycle_count"]/transducer["center_frequency"]) # Pulse length (seconds) given # of cycles
    pulseSteps = int(np.ceil(pulseLen/dt))
    phase_rads = transducer["phase"]*np.pi/180 

    wf_tVec = np.arange(0, (pulseSteps)*dt, dt)
    wf = (np.sin(wf_tVec*transducer["center_frequency"]*2*np.pi+phase_rads)).astype(np.double) 
    window = signal.hann(pulseSteps)
    wf = np.multiply(wf, window) 
    wf = wf/wf.max()

    array = numpy_to_vtk(np.array([1]), deep=True)
    name = "NTransducers"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(np.array([transducer["location"][0]]), deep=True)
    name = "XD0/XCenter"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(np.array([transducer["location"][1]]), deep=True)
    name = "XD0/YCenter"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(np.array([transducer["location"][2]]), deep=True)
    name = "XD0/ZCenter"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(np.array([transducer["radius"]]), deep=True)
    name = "XD0/Radius"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(wf_tVec, deep=True)
    name = "XD0/SignalTimes"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(wf, deep=True)
    name = "XD0/SignalValues"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    ts = tf
    tf = time.time()
    print("Data Fields made in %f seconds!" % (tf-ts))

    ## Writing to File
    writer = vtk.vtkXMLDataSetWriter()
    writer.SetFileName("Composite_Laminate_Example.vtu")
    writer.SetInputData(ugrid)
    writer.Write()
    ts = tf
    tf = time.time()
    print("File Written in %f seconds!" % (tf-ts))

    return

def main():
    #dimentions in m, s, Pa, kg, degrees, except write_interval which is in steps
    lamina_stiffness_matrix = np.array([[174.89e9,     4.09e9,     5.01e9,     0,          0,          0],
                                        [4.09e9,    15.03e9,    5.01e9,     0,          0,          0],
                                        [5.01e9,    5.01e9,     15.03e9,    0,          0,          0],
                                        [0,         0,          0,          5.01e9,     0,          0],
                                        [0,         0,          0,          0,          6.26e9,     0],
                                        [0,         0,          0,          0,          0,          6.26e9]])

    laminate_geometry = {
        "lamiant_dimensions": np.array([180e-3, 180e-3]), 
        "cell_dimensions": np.array([0.36e-3, 0.36e-3, 0.06e-3]), 
        "ply_thickness": 0.12e-3, 
        "ply_orientations": np.array([0, 45, -45, 90, 90, -45, 45, 0])}

    material_properties = {
        "density": 1570, 
        "lamina_stiffness_matrix": lamina_stiffness_matrix}

    timing = {
        "dt": 10e-9, 
        "write_interval": 50, 
        "simulation_duration": 1e-3}

    transducer =  {
        "radius": 2e-3, 
        "center_frequency": 200e3, 
        "cycle_count": 5, 
        "phase": 0, 
        "location": np.array([90e-3, 90e-3, 0.96e-3])}
    
    makeLaminate(laminate_geometry, material_properties, timing, transducer)
    return


if __name__ == '__main__':
    main()