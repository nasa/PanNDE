
import sys
import time
import numpy as np
import process


data_path = sys.argv[1]
data_root = data_path.split('/')[-1]
meta_data_file = sys.argv[2]

print(data_path)
print(data_root)

ts1 = time.time()
ts = ts1
#data_names = ["Sxx","Syy","Szz","Sxy","Syz", "Sxz",]
data_names = ["Sxx","Syy","Szz","Sxy","Syz", "Sxz",]


tStart = 0
time_step_multiplier = 1
loaded_data=np.load(meta_data_file)
tVec = loaded_data["tVec"]
Nt = tVec.shape[0]
ds = loaded_data["ds"]

slice_normal = 2 # 0=X, 1=Y, 2=Z
plane_origin = ds[slice_normal]/100 # offset from the 0 slice

dataset = []
rc = process.slice(data_path, slice_normal, plane_origin=plane_origin)

for data_name in data_names:
	print(data_name)
	frames = rc.createTimeSeriesOfSlice(Nt, tStart, time_step_multiplier, data_name)
	dataset.append(frames)
	tf = time.time()
	print("Data %s completed in %f seconds!" % (data_name,tf-ts))
	ts = tf
	pass

data_file = "%s.npz" % data_root
data_set = np.array(dataset)
np.savez(data_file, data=dataset, data_names=data_names)

tf = time.time()
print("All done in %f seconds!" % (tf-ts1))

