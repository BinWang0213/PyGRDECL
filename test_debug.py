from GRDECL2VTK import *

# 1.1 Set physical dimensions in physDims (m)
physDims=[2000.0,1000,500]

# 1.2 Set model grid dimensions in gridDims
Nx=20;  Ny=20;  Nz=Nx
gridDims=[Nx,Ny,Nz]

# 1.3 Set corner point grid options
opt=buildCPGGrid_opt(disturbed=False,     flat=True,\
                     fault_drop=400 , fault_nx=Nx//2)

# 1.4 Create empty GeologyModel - Build CPG
Model=GeologyModel()
Model.buildCPGGrid(physDims,gridDims,opt)

# 1.5 Compute First TPFA (block centered) Pressure values
# Model.compute_TPFA_Pressure(Press_inj=1,direction="i",Fault_opt=opt )
# Model.plot_scalar("PORO",ITK=True).show()


from utils import *

# 2.1 Create random perm field with normal distribution for each layer
K_LayerPerm=[100,1000,10]
# K,phi=logNormLayers(gridDims,K_LayerPerm,poro_const=0.05)
K,phi=logNormLayers_basc(gridDims,K_LayerPerm,poro_const=0.05)

# 2.2 Update porosity/permeability fields
Update_fields=["PORO","PERMX","PERMY","PERMZ"]
Update_values=[ phi  , K     , K     , 0.1*K ]
Model.UpdateListCellData(var_list=Update_fields,array_list=Update_values)

# 2.3 Compute TPFA (block centered) Pressure values
# Model.compute_TPFA_Pressure(Press_inj=1,direction="i",Fault_opt=opt )
# Model.plot_scalar("Pressure",ITK=True).show()


# 3.1 Set coarsening factor (grid dimensions of coarse cells)
Model.GRDECL_Data.coarse2fine_ratio=[3]*3

# 3.2 Create coarse grid and upscale porosity
Model2=Model.create_coarse_model()
# Model.Upscale_Perm('TPFA_loc')
# Model2.plot_scalar("PORO",ITK=True).show(True)

# 3.3 Upscaling
Model.Upscaler.nlayer=1
Model.Upscale_Perm('TPFA_loc')

# 3.3 Compute Pressure for coarse model
Model2.compute_TPFA_Pressure(Press_inj=1,direction="i",Fault_opt=opt )
# # Model.two_plots_scalar("Pressure").show(False)
# Model2.plot_scalar("Pressure").show()

