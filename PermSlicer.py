import numpy as np
from GRDECL2VTK import * #GRDECL_Reader,GRID2VTU,getIJK


#Import Mesh from SeniorDesign#Import 
Grid1=GRDECL_Viewer(filename='./Example/SPE10B.GRDECL',nx=60,ny=220,nz=85)
#Choose the subset you want to cut
SubsetGrid=Grid1.field_cutter(nx_range=(1,30),ny_range=(1,110),nz_range=(1,40))
#Output to ASCii fiel
SubsetGrid.write_Ascii("./Example/SPE10B")
#Output to Paraview
SubsetGrid.write_VTU(filename='./Example/SPE10B_PermPoro_Subset')
