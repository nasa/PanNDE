import sys
import os
import numpy as np 
import vtk
from vtk.util.numpy_support import numpy_to_vtk
import xml.etree.ElementTree as etree


def get_pvtt_data(file):
	"""
    gathers all the data from a .pvtt file

    Args:
        file: .pvtt file 
    Returns:
        list of all the numpy arrays for every subfile. Every subfile will have the same number of columns but not the same number of rows.
    """
	tree = etree.parse(file)
	pieces = tree.findall('PTable/Piece')
	sub_files = []
	for piece in pieces:
		sub_files.append(piece.attrib['Source'])
	table = []
	for i, file in enumerate(sub_files):
		tblsread=vtk.vtkXMLTableReader()
		tblsread.SetFileName(file)
		tblsread.Update()
		tbl=tblsread.GetOutput()
		N=tbl.GetNumberOfRows()
		M=tbl.GetNumberOfColumns()
		table.append(np.zeros((N,M)))
		for j in range(N):
		    d=tbl.GetRow(j);
		    for k in range(M):
		    	table[i][j,k] = d.GetValue(k).ToFloat()
	table = np.concatenate(table, axis=0)
	return table 



def main(argv):
	"""
	Searches the given directory and its subdirectories for files with a 
		.pvtt extension. Saves all the tables to a .npz file and if the pvtt files includes a none_coords.pvtt file it will eliminates duplicate nodes. 

	Args:
		argv[0]: output file name

	Returns:
		nothing 
	"""

	output_filename = argv[0]
	if len(argv)> 1:
		directory = argv[1]
	else:
		directory = './'

	node_coords = 'node_coords.pvtt'
	node_coords_test = False
	saving_dict={}

	pvtt_files = []
	for file in os.listdir(directory):
		print(file)
		if file.split('.')[-1] == 'pvtt':
			if file == node_coords:
				node_coords_test = True 
			else:
				pvtt_files.append(file)
	if node_coords_test:             
		coords_table = get_pvtt_data(os.path.join(directory, node_coords))
		uniques = []
		for i in range(coords_table.shape[0]):
			dists = []
			for j in range(i):
				dists.append((coords_table[i][0] - coords_table[j][0])**2+(coords_table[i][1] - coords_table[j][1])**2+ (coords_table[i][2] - coords_table[j][2])**2)
			if all(np.array(dists)>0.0):
				uniques.append(i)
		saving_dict['node_coordinates'] = coords_table[uniques,:]
	else:
		print("Cannot find %s in %s" % (node_coords, directory))
		pass


	variables = []
	timesteps = []
	tables=[]
	for file in pvtt_files:
		parts = file.split('.')[0].split('_')
		if len(parts)>2:
			variable = '_'.join(parts[:-1])
		else:
			variable = parts[0]
		timestep = int(parts[-1])
		if variable not in variables:
			variables.append(variable)
			timesteps.append([])
			tables.append([]) 
		idx = variables.index(variable)
		timesteps[idx].append(timestep)
		tables[idx].append(get_pvtt_data(os.path.join(directory, file)))

	for i, variable in enumerate(variables):
		var_ts = np.array(timesteps[i])
		num_ts = var_ts.max()+1
		num_elements = tables[i][0].shape[0]
		if node_coords_test:
			num_elements = len(uniques)
		num_variable_components = tables[i][0].shape[1]

		sorted_data = np.zeros((num_elements, num_ts, num_variable_components))
		for t in var_ts:
			if node_coords_test:
				sorted_data[:,t,:]=tables[i][t][uniques,:]
			else:
				sorted_data[:,t,:]=tables[i][t][:,:]
		saving_dict[variable]=sorted_data

	print("saving to %s" % output_filename)
	np.savez(output_filename, data_dict=saving_dict)
	return True


if __name__ == '__main__':
    main(sys.argv[1:])


### to open:
# file = np.load(output_filename, allow_pickle=True)
# data_dict = file['data_dict'].item()
