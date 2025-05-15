import vtk
import numpy as np


def saveSlice(fname_base, timestep, data_sets, resample_sz, img_dims, sliceNormal, sliceOrigin)
:
	timestep = i+1+tstart
	filename = "%s%i.pvtu"%(fname_base, timestep)
	print(filename)
	reader = vtk.vtkXMLPUnstructuredGridReader()
	reader.SetFileName(filename)
	reader.Update()

	pt2cell = vtk.vtkPointDataToCellData()
	pt2cell.SetInputConnection(reader.GetOutputPort())
	pt2cell.Update()

	plane = vtk.vtkPlane()
	plane.SetNormal(sliceNormal[0], sliceNormal[1], sliceNormal[2])
	plane.SetOrigin(sliceOrigin[0], sliceOrigin[1], sliceOrigin[2])
  
	slicer = vtk.vtkCutter()#vtk.vtkPlaneCutter()
	slicer.SetInputConnection(pt2cell.GetOutputPort())
	slicer.SetCutFunction(plane)
	slicer.Update()

	pdat = slicer.GetOutput()

	resampler = vtk.vtkResampleToImage()
	resampler.SetInputDataObject(pdat)
	resampler.SetSamplingDimensions(int(resample_sz[0]), int(resample_sz[1]), int(resample_sz[2]))
	resampler.Update()
	img_data = resampler.GetOutput()

	for data_name in data_sets:
		frame = np.zeros((img_dims[0], img_dims[1])).T
		Result_raw = img_data.GetPointData().GetAbstractArray(data_name)
		for kx in range(img_dims[0]): 
			for ky in range(img_dims[1]): 
				frame[ky, kx] = Result_raw.GetTuple1(kx+ky*img_dims[0])
				pass
			pass
		pass
		ofname = '%s_%04d.npz' %(data_name, timestep)
		print(ofname)
		np.savez(ofname, data=frame)
		pass
	return 1

def saveSlice0(fname_base, timestep, data_sets, sliceNormal, sliceOrigin):
	filename = "%s%i.pvtu"%(fname_base, timestep)
	print(filename)

	reader = vtk.vtkXMLPUnstructuredGridReader()
	reader.SetFileName(filename)
	reader.Update()

	pt2cell = vtk.vtkPointDataToCellData()
	pt2cell.SetInputConnection(reader.GetOutputPort())
	pt2cell.Update()

	plane = vtk.vtkPlane()
	plane.SetNormal(sliceNormal[0], sliceNormal[1], sliceNormal[2])
	plane.SetOrigin(sliceOrigin[0], sliceOrigin[1], sliceOrigin[2])

	slicer = vtk.vtkCutter()#vtk.vtkPlaneCutter()
	slicer.SetInputConnection(pt2cell.GetOutputPort())
	slicer.SetCutFunction(plane)
	slicer.Update()

	pdat = slicer.GetOutput()
	
	bounds = pdat.GetBounds()
	dims = np.array([bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4]])
	cell_aabb = pdat.GetCell(5).GetBounds()

	resample_sz = np.zeros(3)

	resample_sz[0] = int(dims[0]/(cell_aabb[1]-cell_aabb[0]))
	resample_sz[1]= int(dims[1]/(cell_aabb[3]-cell_aabb[2])) 
	resample_sz[2] = 1

	resampler = vtk.vtkResampleToImage()
	resampler.SetInputDataObject(pdat)
	resampler.SetSamplingDimensions(int(resample_sz[0]), int(resample_sz[1]), int(resample_sz[2]))
	resampler.Update()
	img_data = resampler.GetOutput()
	img_dims = img_data.GetDimensions()
	#spacing = img_data.GetSpacing()

	for data_name in data_sets:
		frame = np.zeros((img_dims[0], img_dims[1])).T
		Result_raw = img_data.GetPointData().GetAbstractArray(data_name)
		for kx in range(img_dims[0]): 
			for ky in range(img_dims[1]): 
				frame[ky, kx] = Result_raw.GetTuple1(kx+ky*img_dims[0])
				pass
			pass
		pass
		ofname = '%s_%04d.npz' %(data_name, timestep)
		print(ofname)
		np.savez(ofname, data=frame)
		pass
	return resample_sz, img_dims



fname_base = 'Putput_file_base_'
data_sets = ['Vx', 'Vy', 'Vz']
tstart = 0
Nt = 11

sliceNormal = [0,0,1] # [0,0,1] = Z-normal
sliceOrigin = [0,0,1.0e-5] # where in space is the plane's origin 

timestep = tstart
resample_sz, img_dims = saveSlice0(fname_base, timestep, data_sets, sliceNormal, sliceOrigin)

for i in range(Nt-1):
	timestep = i+1+tstart
	saveSlice(fname_base, timestep, data_sets, resample_sz, img_dims, sliceNormal, sliceOrigin)

	pass

