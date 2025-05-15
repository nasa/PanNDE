#!#!/usr/bin/env python

import sys
import time
import math
import itertools
import numpy as np
import scipy.io as sio
from scipy import signal
import vtk
from vtk.util.numpy_support import numpy_to_vtk


def makePlate(ds, len_s, hole_locs, hole_rads):
    """
    Create rectangular prism geometry

    Parameters
    ----------
    ds      : a three element array of floats
        the cell dimensions in meters
    len_s   : a three element array of floats
        the encompassing volume dimensions in meters


    Returns
    -------
    Geom    : a 3D numpy array of 0 zeros and ones
        The Geometry of the dog bone inside the array defined by 1s and 0s are empty space

    """

    # Solid cells
    Ng = [int(np.rint(len_s[a]/ds[a])) for a in range(3)]

    Geom = np.zeros((Ng[0]+2, Ng[1]+2, Ng[2]+2), dtype=np.uint8)
    Geom[1:-1, 1:-1, 1:-1] = 1

    print("Geom shape: %dx%dx%d" %(Ng[0]+2, Ng[1]+2, Ng[2]+2))


    for hole, hole_loc in enumerate(hole_locs):
        rad = hole_rads[hole]
        print("placing hole %d" %(hole+1))
        loc_N = [int(np.rint(hole_loc[a]/ds[a])) for a in range(2)]
        print("Hole located at %d, %d" %(loc_N[0], loc_N[1]))
        rad_Ng = int(np.rint(rad/ds[0]))
        print("radius = %d" % (rad_Ng))
        hole_cube = np.zeros((2*rad_Ng, 2*rad_Ng, Ng[2]+2), dtype=np.uint8)
        hole_cube[:,:, 1:-1]=1
        for i in range(hole_cube.shape[0]):
            for j in range(hole_cube.shape[1]):
                dist = ((i-rad_Ng+.5))**2+((j-rad_Ng+.5))**2
                if dist < rad_Ng**2:
                    hole_cube[i,j,:] = 0
                    pass 
                pass 
            pass 
        Geom[loc_N[0]-rad_Ng+1:loc_N[0]+rad_Ng+1,loc_N[1]-rad_Ng+1:loc_N[1]+rad_Ng+1, : ] = hole_cube
        pass 

    return Geom


# Plate
def main(argv):
    """
    main, runs the script

    Parameters
    ----------
    NONE

    Returns
    -------
    NONE
    """

    # flag = argv[0]
    # file_root = argv[1]
    # metadatFilename = argv[2]
    # num_sims = int(argv[3])

    #plate_dims = [416.56e-3, 144.0e-3, 0.81e-3] #m
    plate_dims = [208e-3, 37.0e-3, 0.81e-3] #m
    ds = [2.5e-4, 2.5e-4, 1.0e-5] #m

    x_cut = 50e-3

    y_offset = 24.24e-3
    x_offset = 103.75e-3 - x_cut #128.65e-3 - 9.65e-3 - 23.88e-3 
    y_spacing = 23.88e-3
    x_spacing = 15.25e-3

    rows = 1
    cols = 2

    hole_locs = []
    for i in range(cols):
        for j in range(rows):
            hole_locs.append((x_offset+x_spacing*i, y_offset+y_spacing*j))
            print("Hole %d located at %f, %f mm" %(i*rows+j+1, hole_locs[-1][0],hole_locs[-1][1] ))
            pass 
        pass 

    hole_rads = [4.0386e-3 for i in range(rows*cols)]


    tStart = 0.0
    dt = (0.125/12)*1e-7
    dt_out = (0.125/3)*1e-6 # s
    num_outputs = 6
    tStop = dt_out*num_outputs
    #tStop = 166.625e-6

    tVec = np.arange(tStart, tStop, dt)
    wo_arr = np.arange(tStart, tStop, dt_out)
    fname = 'small_AL_plate_w_holes_Transducer_Example.vtu'
    flag= makePlateSim(fname, plate_dims, ds, hole_locs, hole_rads, tVec, wo_arr)
    return 


def smoothPulse(t,pulseWidth=1.0):
  return 0.5*(1.0-np.cos(2.0*math.pi*t/pulseWidth))*(t<pulseWidth)

def makePlateSim(fname, plate_dim, ds, hole_locs, hole_rads, tVec, wo_arr):
    """
    Creates the VTU File for a specified geometry 

    Parameters
    ----------
    fname      : a string
        the vtu file name
    block_dim  : an array of float
        the dimensions of the block in meters
    ds         : an array of floats with the dimensions of each cell
    vol_frac   : target porosity volume fraction

    Returns
    -------
    NONE
    """
    ts1 = time.time()
    ts = ts1


    dt = tVec[1]-tVec[0]

    x_cut = 50e-3

    xdcr_x = 200.275e-3 - x_cut # 
    xdcr_y = 18e-3 # 
    xdcr_z = ds[2]*0.01
    xdcr_r = 6.35e-3

    # material(s) definitions
    # 2024-T4
    num_mat = 1
    rho = [2780] # kg/m^3
    E = [73.1e9] # Pa
    G = [28.9e9] # Pa
    nu = [0.33]
    C = np.zeros((21, num_mat))
    for i in range(num_mat):
        C[0, i] = E[i]*(1-nu[i])/((nu[i]+1)*(1-2*nu[i]))
        C[1, i] = E[i]*nu[i]/((nu[i]+1)*(1-2*nu[i]))
        C[2, i] = E[i]*nu[i]/((nu[i]+1)*(1-2*nu[i]))
        C[6, i] = E[i]*(1-nu[i])/((nu[i]+1)*(1-2*nu[i]))
        C[7, i] = E[i]*nu[i]/((nu[i]+1)*(1-2*nu[i]))
        C[11, i] = E[i]*(1-nu[i])/((nu[i]+1)*(1-2*nu[i]))
        C[15, i] = G[i]
        C[18, i] = G[i]
        C[20, i] = G[i]
        pass

    len_s = plate_dim
    Ng = [int(np.rint(len_s[a]/ds[a])) for a in range(3)]
    Geom = makePlate(ds, len_s, hole_locs, hole_rads)
    tf = time.time()
    print("Geometry Made in %f seconds!" % (tf-ts))

    numCells = np.sum(Geom)
    points = vtk.vtkPoints()
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.Allocate(numCells)
    pts = np.ones((Ng[0]+1, Ng[1]+1, Ng[2]+1))*-1
    rho_arr = np.ones(numCells)
    C_arr = np.ones((21, numCells))

    ts = tf
    tf = time.time()
    print("arrays Made in %f seconds!" % (tf-ts))


    pindx = 0
    cell_idx = 0
    for i, j, k in itertools.product(range(Ng[0]+1), range(Ng[1]+1), range(Ng[2]+1)):
        test = np.sum(Geom[i:i+2, j:j+2, k:k+2])
        if test >= 1:
            pts[i, j, k] = pindx
            pindx = pindx + 1
            points.InsertPoint(int(pts[i, j, k]), [i*ds[0], j*ds[1], k*ds[2]])
            pass
            if Geom[i, j, k] > 0:
                cube = [int(pts[i-1, j-1, k-1]), int(pts[i, j-1, k-1]), int(pts[i, j, k-1]), int(pts[i-1, j, k-1]), int(pts[i-1, j-1, k]), int(pts[i, j-1, k]), int(pts[i, j, k]), int(pts[i-1, j, k])]
                # [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1],[0, 1, 1]]
                if (np.array(cube) < 0).any():
                    print("cube failure at %i,%i,%i\n"%(i, j, k))
                    pass
                ugrid.InsertNextCell(vtk.VTK_HEXAHEDRON, 8, cube)
                C_arr[:, cell_idx] = C[:, Geom[i, j, k]-1]
                rho_arr[cell_idx] = rho[Geom[i, j, k]-1]
                cell_idx = cell_idx + 1
                pass
            pass
        #print("Slice %d of %d done!" % (i, Ng[0]))
        pass

    ts = tf
    tf = time.time()
    print("Unstructured Grid Made in %f seconds!" % (tf-ts))
    ugrid.SetPoints(points)

    C_idx = 0
    for i, j in itertools.combinations_with_replacement(range(6), 2):
        name = "C%d%d" %(i+1, j+1)
        array = numpy_to_vtk(C_arr[C_idx, :], deep=True)
        array.SetName(name)
        ugrid.GetCellData().AddArray(array)
        C_idx = C_idx +1
        pass

    array = numpy_to_vtk(rho_arr, deep=True)
    name = "density"
    array.SetName(name)
    ugrid.GetCellData().AddArray(array)

    array = numpy_to_vtk(np.array([dt]), deep=True)
    name = "dt"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(wo_arr, deep=True)
    name = "write_times"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

#### Transducers
    array = numpy_to_vtk(np.array([3]), deep=True)
    name = "NTransducers"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)
# Driving Transducer
 # Windowed 
    domFreq = 250e3  # Hz
    numCycles = 5
    pulseLen = (numCycles/domFreq) # Pulse length (seconds) given # of cycles
    pulseSteps = int(np.ceil(pulseLen/dt)) # # of pulse time steps 

    amp = 1
    wf_tVec = np.arange(0, (pulseSteps)*dt, dt)
    wf = (amp*np.sin(wf_tVec*domFreq*2*np.pi)).astype(np.double) 
    window = signal.windows.hann(pulseSteps)
    wf = np.multiply(wf, window) 
    wf = wf/wf.max()

    array = numpy_to_vtk(np.array([xdcr_x]), deep=True)
    name = "XD0/XCenter"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(np.array([xdcr_y]), deep=True)
    name = "XD0/YCenter"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(np.array([xdcr_z]), deep=True)
    name = "XD0/ZCenter"
    array.SetName(name)
    ugrid.GetFieldData().AddArray(array)

    array = numpy_to_vtk(np.array([xdcr_r]), deep=True)
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

    # Recieving Transducer 

    xdcr_r = 3.175e-3
    xdcr_x0 = 111.375e-3 - x_cut # 
    xdcr_x_spacing = 25e-3
    xdcr_y0 = 14e-3 # 
    xdcr_y_spacing = 24e-3
    xdcr_z = ds[2]*0.01

    wf_tVec_pass = wf_tVec[:10]
    wf_pass = np.zeros(wf_tVec_pass.shape)

    i = 1
    for a in range(2):
        for b in range(1):
            xdcr_x = xdcr_x0+a*xdcr_x_spacing
            xdcr_y = xdcr_y0+b*xdcr_y_spacing

            array = numpy_to_vtk(np.array([xdcr_x]), deep=True)
            name = "XD%d/XCenter" % i
            array.SetName(name)
            ugrid.GetFieldData().AddArray(array)

            array = numpy_to_vtk(np.array([xdcr_y]), deep=True)
            name = "XD%d/YCenter" % i
            array.SetName(name)
            ugrid.GetFieldData().AddArray(array)

            array = numpy_to_vtk(np.array([xdcr_z]), deep=True)
            name = "XD%d/ZCenter" % i
            array.SetName(name)
            ugrid.GetFieldData().AddArray(array)

            array = numpy_to_vtk(np.array([xdcr_r]), deep=True)
            name = "XD%d/Radius" % i
            array.SetName(name)
            ugrid.GetFieldData().AddArray(array)

            array = numpy_to_vtk(wf_tVec_pass, deep=True)
            name = "XD%d/SignalTimes" % i
            array.SetName(name)
            ugrid.GetFieldData().AddArray(array)

            array = numpy_to_vtk(wf_pass, deep=True)
            name = "XD%d/SignalValues" % i
            array.SetName(name)
            ugrid.GetFieldData().AddArray(array)
            i = i+1
            pass 
        pass 


    
    ts = tf
    tf = time.time()
    print("Data Fields made in %f seconds!" % (tf-ts))

    writer = vtk.vtkXMLDataSetWriter()
    writer.SetFileName(fname)
    writer.SetInputData(ugrid)
    writer.Write()
    ts = tf
    tf = time.time()
    print("File Written in %f seconds!" % (tf-ts))
    print("Whole Script Executed in %f seconds!" % (tf-ts1))

    print("Ng = [%d, %d, %d] with ds = [%f, %f, %f] mm" % (Ng[0], Ng[1], Ng[2], ds[0]*1000, ds[1]*1000, ds[2]*1000))
    print("Total Cells = %d" % numCells)
    print("Total Cells Created = %d" % ugrid.GetNumberOfCells())
    print("Total Points = %d" % ugrid.GetNumberOfPoints())

    print("Real ds dimension is [%f, %f, %f] m with an error of [%f, %f, %f] m" % (Ng[0]*ds[0], Ng[1]*ds[1], Ng[2]*ds[2], len_s[0]-Ng[0]*ds[0], len_s[1]-Ng[1]*ds[1], len_s[2]-Ng[2]*ds[2]))
    print("Written to %s." % fname)
    return 1

if __name__ == '__main__':
    main(sys.argv[1:])
