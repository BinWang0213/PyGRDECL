from utils import *

class Upscaler:
    def __init__(self):
        """
        Fields and methods to upscale a geomodel
        """
        self.Coarse_Mod=None
        self.Coarse_Mod_Partition=[]
        self.local_Mod=None

    def Upscale_Arithmetic_mean(self,FineModel, scalar):
        Finefield = FineModel.GRDECL_Data.SpatialDatas[scalar]
        FGrid = FineModel.GRDECL_Data
        CGrid = self.Coarse_Mod.GRDECL_Data
        vol = FGrid.VOL
        P0 = self.Coarse_Mod_Partition
        # 1. Arithmetic Average (Weighted by volumes) per coarse grid cells
        CoarseField = np.bincount(P0, weights=vol * Finefield) / np.bincount(P0, weights=vol)
        self.Coarse_Mod.UpdateCellData(varname=scalar, array=CoarseField)
        m = np.min(CoarseField)
        M = np.max(CoarseField)
        plot_hist(self.Coarse_Mod.GRDECL_Data.SpatialDatas[scalar], \
                  varname=scalar + ": " + "arithmetic: min %3f - max %3f" % (m, M))
