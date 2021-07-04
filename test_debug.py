from GRDECL2VTK import *
import numpy as np, math

FilePath="ExampleData/"
Basefilename="CPGrid00"
FileExtension="GRDECL"

Model=GeologyModel()
Model.fname=FilePath+Basefilename + "." + FileExtension

Nxyz=Nx=Ny=Nz=4
divisor=2
xmax=1
ymax=xmax
zmax=xmax/divisor
Nz//=divisor
Nz=3
opt={'disturbed':False,'flat':True}
Model.buildCPGGrid([xmax,ymax,zmax],[Nx,Ny,Nz],opt,faultDrop=0.15)
Model.processGRDECL()
