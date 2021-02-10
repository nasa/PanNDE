# Python based vtk data reconstructor with frequency domain post processing. 
#  - included for the V&V work performed in the paper, and as an example  
#    post processing utility
#
# Notices:
# Copyright 2021 United States Government as represented by the 
# Administrator of the National Aeronautics and Space Administration. 
# No copyright is claimed in the United States under Title 17, U.S. Code. 
# All Other Rights Reserved.
# 
# Googletest is a product of Google Inc. and is subject to the following:
# 
# Copyright 2008, Google Inc. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met:
#  * Redistributions of source code must retain the above copyright notice, 
#    this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice, 
#    this list of conditions and the following disclaimer in the documentation 
#    and/or other materials provided with the distribution.
#  * Neither the name of Google Inc. nor the names of its contributors may be used 
#    to endorse or promote products derived from this software without specific 
#    prior written permission.
#
# Disclaimers
# No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, 
# EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT 
# THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT 
# SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO 
# THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY 
# GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE 
# PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT 
# AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE 
# ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS." 
#
# Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES 
# GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S 
# USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES 
# ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S 
# USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, 
# ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  
# RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS 
# AGREEMENT.


import vtk
import scipy.fft as spfft
import scipy.signal as spsig
import numpy as np

class reconstructor:
  def __init__(self,fname_base=""):
    self.setFilenameBase(fname_base)
    return

  def setFilenameBase(self,fname_base):
    self._fname_root=fname_base
    return

  def generate3DFFT(self,Nt,tstart,dt,z=0.,data_name="Vz",time_step_multiplier=1):
    tukey_time=spsig.tukey(Nt,alpha=0.5)
    for kt in range(Nt):
      tidx=(tstart+kt)*time_step_multiplier
      data,dims,ds=self._loadTimeStepDataSlice(tidx,z,data_name=data_name)
      if kt==0:
        result=np.zeros([data.shape[0],data.shape[1],Nt],dtype=np.cdouble)
        tukey_space=self._makeTukeySpatialWindow(dims,alpha=0.5)
      result[:,:,kt]=data*tukey_space*tukey_time[kt]
    wavenum_x,wavenum_y=self._makeWavenumbers(dims,ds)
    frequency=spfft.fftshift(spfft.fftfreq(Nt,dt))
    result=spfft.fftshift(spfft.fftn(result))
    return result,wavenum_x,wavenum_y,frequency

  def generateWavenumberDataAtFrequency(self,Nt,tstart,freq,dt,z=0.,
                                        data_name="Vz",time_step_multiplier=1,xdrad=0.002):
    tukey_time=spsig.tukey(Nt,alpha=0.5)
    for kt in range(Nt):
      tidx=(tstart+kt)*time_step_multiplier
      data,dims,ds=self._loadTimeStepDataSlice(tidx,z,data_name=data_name)
      if kt==0:
        Acc=np.zeros(data.shape,dtype=np.cdouble)
      Acc+=data*np.exp(-1j*2*np.pi*freq*tidx*dt)*dt*tukey_time[kt]
    wavenum_x,wavenum_y=self._makeWavenumbers(dims,ds)
    #tukey_space=self._makeTukeySpatialWindow(dims,alpha=0.5)
    tukey_space=self._makeNotchedTukeySpatialWindow(dims,ds,xdrad=xdrad,alpha=0.5)
    
    result=spfft.fftshift(spfft.fftn(Acc*tukey_space))
    return result,wavenum_x,wavenum_y,Acc

  # def writeSegment(self,filename,segment,timestep,z=0.):
  #   dat=self._readSegmentAtTimestep(segment,timestep,z)
  #   writer=vtk.vtkXMLImageDataWriter()
  #   writer.SetFileName(filename)
  #   writer.SetInputData(dat)
  #   writer.Update()
  #   return
  
  # def writeSlices(self,Nt,write_root,read_root,z=0.):
  #   for kt in range(Nt):
  #     wname="%s_%i.vti"%(write_root,kt)
  #     rname="%s_%i.pvtu"%(read_root,kt)
  #     self.writeSlice(wname,rname,z)
  #   return

  # def writeSlice(self,filename,readname,z=0.):
  #   dat=self._readSliceFromVTUFile(readname,z)
  #   img=self._makeVTKImageFromPolydata(dat)
  #   writer=vtk.vtkXMLImageDataWriter()
  #   writer.SetFileName(filename)
  #   writer.SetInputData(img)
  #   writer.Write()
  #   return

## private stuff
  def _makeWavenumbers(self,dims,spacing):
    wavenum_x=spfft.fftshift(spfft.fftfreq(dims[0],spacing[0]))
    wavenum_y=spfft.fftshift(spfft.fftfreq(dims[1],spacing[1]))
    return wavenum_x,wavenum_y

  def _makeTukeySpatialWindow(self,dims,alpha=0.5):
    tukey_win_x=spsig.tukey(dims[0],alpha=alpha)
    tukey_win_y=spsig.tukey(dims[1],alpha=alpha)
    tukey_space=np.outer(tukey_win_y,tukey_win_x)
    return tukey_space

  def _makeNotchedTukeySpatialWindow(self,dims,ds,xdrad=0.002,alpha=0.5):
    tukey_win_x=spsig.tukey(dims[0],alpha=alpha)
    tukey_win_y=spsig.tukey(dims[1],alpha=alpha)

    x=np.linspace(0,dims[0]-1,dims[0])*ds[0]
    y=np.linspace(0,dims[1]-1,dims[1])*ds[1]

    xc=x[-1]/2.
    x=x-xc
    yc=y[-1]/2.
    y=y-yc
    X,Y=np.meshgrid(x,y)
    R=np.sqrt(X*X+Y*Y)
    notch_win=((R<xdrad)+0.5*(R>xdrad)*(1.+np.cos(np.pi*(R-xdrad)/xdrad)))*(R>=0.)*(R<2.*xdrad)
    notch_win=1.*(R<xdrad)+0.5*(R>=xdrad)*(1.+np.cos(np.pi*(R-xdrad)/xdrad))*(R<(2.*xdrad))

    tukey_space=np.outer(tukey_win_y,tukey_win_x)-notch_win
    return tukey_space

  def _loadTimeStepDataSlice(self,tidx,z,data_name="Vz"):
    print("load timestep: %i"%tidx)
    pdat=self._readSegmentAtTimestep(tidx,z)
    img_data=self._makeVTKImageFromPolydata(pdat)
    data=self._makeNumpySliceFromImg(img_data,data_name=data_name)
    dims=img_data.GetDimensions()
    spacing=img_data.GetSpacing()
    return data,dims,spacing

  def _makeNumpySliceFromImg(self,img_data,data_name="Vz"):
    dims=img_data.GetDimensions()
    Result=np.zeros(dims[0:2]).T
    Result_raw=img_data.GetPointData().GetAbstractArray(data_name)
    for kx in range(dims[0]): 
      for ky in range(dims[1]): 
        Result[ky,kx]=Result_raw.GetTuple1(kx+ky*dims[0])
    return Result

  def _readSliceFromVTUFile(self,filename,z=0.):
    # Read the source file.
    reader=vtk.vtkXMLPUnstructuredGridReader()#vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()  # Needed because of GetScalarRange    
    pt2cell=vtk.vtkPointDataToCellData()
    pt2cell.SetInputConnection(reader.GetOutputPort())
    pt2cell.Update()
    plane=vtk.vtkPlane()
    plane.SetNormal(0.,0.,1.)
    plane.SetOrigin(0.,0.,z)
    slicer=vtk.vtkCutter()#vtk.vtkPlaneCutter()
    slicer.SetInputConnection(pt2cell.GetOutputPort())
    slicer.SetCutFunction(plane)#SetPlane(plane)#
    slicer.Update()
    output=slicer.GetOutput()
    return output

  def _makeVTKImageFromPolydata(self,dat):
    bounds=dat.GetBounds()
    dims=np.array([bounds[1]-bounds[0],bounds[3]-bounds[2]])
    maxdim=np.max(dims)
    cell_aabb=dat.GetCell(0).GetBounds()
    dx=cell_aabb[1]-cell_aabb[0]
    dy=cell_aabb[3]-cell_aabb[2]
    resampler=vtk.vtkResampleToImage()
    resampler.SetInputDataObject(dat)
    resampler.SetSamplingDimensions(int(dims[0]/dx),int(dims[1]/dy),1)
    resampler.Update()
    output=resampler.GetOutput()
    return output

  def _readSegmentAtTimestep(self,timestep,z=0.):
    filename="%s_%i.pvtu"%(self._fname_root,timestep)
    output=self._readSliceFromVTUFile(filename,z)
    return output
