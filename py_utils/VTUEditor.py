import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk


def get_field_data_array_names(filename):
	reader = vtk.vtkXMLUnstructuredGridReader()
	reader.SetFileName(filename)
	reader.Update()
	ugrid = reader.GetOutput()
	fields = ugrid.GetFieldData()
	num_arrays = fields.GetNumberOfArrays()
	field_data_names = []
	for i in range(num_arrays):
		field_data_names.append(fields.GetArray(i).GetName())
		pass
	return field_data_names

def get_cell_data_array_names(filename):
	reader = vtk.vtkXMLUnstructuredGridReader()
	reader.SetFileName(filename)
	reader.Update()
	ugrid = reader.GetOutput()
	cells = ugrid.GetCellData()
	num_arrays = cells.GetNumberOfArrays()
	cell_data_names = []
	for i in range(num_arrays):
		cell_data_names.append(cells.GetArray(i).GetName())
		pass
	return cell_data_names

def edit_field_data_array(in_filename, out_filename, field_data_array_to_replace, new_array):
	reader = vtk.vtkXMLUnstructuredGridReader()
	reader.SetFileName(in_filename)
	reader.Update()
	ugrid = reader.GetOutput()
	fields = ugrid.GetFieldData()
	if 1==fields.HasArray(field_data_array_to_replace):  
		fields.RemoveArray(field_data_array_to_replace)
	new_vtk_array=numpy_to_vtk(new_array)
	new_vtk_array.SetName(field_data_array_to_replace)
	fields.AddArray(new_vtk_array)
	ugrid.SetFieldData(fields)
	writer=vtk.vtkXMLUnstructuredGridWriter()
	writer.SetFileName(out_filename)
	writer.Update()
	return

def get_field_data_array_info(in_filename, field_array_data_name):
	reader = vtk.vtkXMLUnstructuredGridReader()
	reader.SetFileName(in_filename)
	reader.Update()
	ugrid = reader.GetOutput()
	meta = ugrid.GetFieldData()
	field_data = meta.GetArray(field_array_data_name)
	print("Field Data Array: %s" % field_array_data_name)
	field_data_type = field_data.GetDataType()
	data_types = ["VTK_TYPE_CHAR_IS_SIGNED or VTK_VOID", "VTK_BIT", "VTK_CHAR", "VTK_UNSIGNED_CHAR", "VTK_SHORT", "VTK_UNSIGNED_SHORT ", 
					"VTK_INT", "VTK_UNSIGNED_INT", "VTK_LONG", "VTK_UNSIGNED_LONG", "VTK_FLOAT", "VTK_DOUBLE", "VTK_ID_TYPE", "VTK_STRING",
					"VTK_OPAQUE", "VTK_SIGNED_CHAR", "VTK_UNSIGNED_LONG_LONG" ]
	if field_data_type > len(data_types)-1:
		print("Data Type: Unknown")
	else:
		print("Data Type: %s" % data_types[field_data_type])
	print("Array Size: %s" %(str(field_data.GetNumberOfTuples())))
	return

def get_cell_data_array_info(in_filename, cell_array_data_name):
	reader = vtk.vtkXMLUnstructuredGridReader()
	reader.SetFileName(in_filename)
	reader.Update()
	ugrid = reader.GetOutput()
	meta = ugrid.GetCellData()
	cell_data = meta.GetArray(cell_array_data_name)
	print("Cell Data Array: %s" % cell_array_data_name)
	cell_data_type = cell_data.GetDataType()
	data_types = ["VTK_TYPE_CHAR_IS_SIGNED or VTK_VOID", "VTK_BIT", "VTK_CHAR", "VTK_UNSIGNED_CHAR", "VTK_SHORT", "VTK_UNSIGNED_SHORT ", 
					"VTK_INT", "VTK_UNSIGNED_INT", "VTK_LONG", "VTK_UNSIGNED_LONG", "VTK_FLOAT", "VTK_DOUBLE", "VTK_ID_TYPE", "VTK_STRING",
					"VTK_OPAQUE", "VTK_SIGNED_CHAR", "VTK_UNSIGNED_LONG_LONG" ]
	if cell_data_type > len(data_types)-1:
		print("Data Type: Unknown")
	else:
		print("Data Type: %s" % data_types[cell_data_type])
	print("Array Size: %s" %(str(cell_data.GetNumberOfTuples())))
	return







