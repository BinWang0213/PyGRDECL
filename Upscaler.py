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

    def Upscale_Geometric_mean(self,FineModel, scalar):
        Finefield=FineModel.GRDECL_Data.SpatialDatas[scalar]
        FGrid=FineModel.GRDECL_Data
        CGrid=self.Coarse_Mod.GRDECL_Data
        vol =FGrid.VOL
        P0=self.Coarse_Mod_Partition
        # 1. Arithmetic Average (Weighted by volumes) per coarse grid cells
        CoarseField = np.exp(np.bincount(P0, weights=vol * np.log(Finefield)) / np.bincount(P0, weights=vol))
        self.Coarse_Mod.UpdateCellData(varname=scalar, array=CoarseField)
        m=np.min(CoarseField)
        M=np.max(CoarseField)
        plot_hist(self.Coarse_Mod.GRDECL_Data.SpatialDatas[scalar],\
                  varname=scalar + ": " + "geometric: min %3f - max %3f" %(m,M))

    def Upscale_Harmonic_mean(self,fineModel, scalar):
        Finefield=fineModel.GRDECL_Data.SpatialDatas[scalar]
        FGrid=fineModel.GRDECL_Data
        CGrid=self.Coarse_Mod.GRDECL_Data
        vol =FGrid.VOL
        P0=self.Coarse_Mod_Partition
        # 2. Harmonic Average (Weighted by volumes) per coarse grid cells
        CoarseField=np.bincount(P0,weights=vol)/np.bincount(P0,weights=vol/Finefield)
        self.Coarse_Mod.UpdateCellData(varname=scalar, array=CoarseField)
        m=np.min(CoarseField)
        M=np.max(CoarseField)
        plot_hist(self.Coarse_Mod.GRDECL_Data.SpatialDatas[scalar],\
                  varname=scalar + ": " + "harmonic: min %3f - max %3f" %(m,M))

    def Upscale_Harmx_mean(self,FineModel, scalar):
        KF=FineModel.GRDECL_Data.SpatialDatas[scalar]
        FGrid = FineModel.GRDECL_Data
        CGrid = self.Coarse_Mod.GRDECL_Data
        vol = FGrid.VOL
        P0 = self.Coarse_Mod_Partition
        # Computing coarse K values
        # In each coarse cell: for each fixed (i,j) compute x harm mean
        # compute arithmetic mean of the harm values
        # 3a. Harmonic(in x)-Arithmetic Average
        KC=np.zeros((CGrid.NX,CGrid.NY,CGrid.NZ),dtype=float)
        Kvol=vol/KF
        for k in range(CGrid.NZ):
            for j in range(CGrid.NY):
                for i in range(CGrid.NX):
                    indices=np.where(P0==i+j*CGrid.NX+k*CGrid.NX*CGrid.NY)
                    # Local & x partition based: KCloc is a Vector of Harm means in local x slices
                    rx = FGrid.Rx[i];
                    ry = FGrid.Ry[j];
                    rz = FGrid.Rz[k];

                    # z only growing local partition vector Plocy
                    Plocx = np.ones(rx * ry * rz, dtype=int)
                    compt = 0
                    for iii in range(0, rx * ry * rz, rx):
                        Plocx[iii:iii + rx] = compt
                        compt += 1

                    vol_loc=np.bincount(Plocx,weights=vol[indices])
                    KCloc=vol_loc/np.bincount(Plocx,weights=Kvol[indices])
                    #Arithm mean of all the local Harm means
                    KC[i,j,k] = np.sum(KCloc*vol_loc)/np.sum(vol_loc)
        self.Coarse_Mod.UpdateCellData(varname=scalar, array=KC.reshape((-1,1),order='F'))
        m=np.min(KC)
        M=np.max(KC)
        plot_hist(self.Coarse_Mod.GRDECL_Data.SpatialDatas[scalar],\
                  varname=scalar + ": " + "harmonicx: min %3f - max %3f" %(m,M))

    def Upscale_Harmy_mean(self,FineModel, scalar):
        KF=FineModel.GRDECL_Data.SpatialDatas[scalar]
        FGrid = FineModel.GRDECL_Data
        CGrid = self.Coarse_Mod.GRDECL_Data
        vol = FGrid.VOL
        P0 = self.Coarse_Mod_Partition
        # Computing coarse K values
        # 3b. Harmonic(in y)-Arithmetic Average
        KC=np.zeros((CGrid.NX,CGrid.NY,CGrid.NZ),dtype=float)
        rx=FGrid.NX//CGrid.NX;
        ry=FGrid.NY//CGrid.NY;
        rz=FGrid.NZ//CGrid.NZ
        Kvol=vol/KF
        for k in range(CGrid.NZ):
            for j in range(CGrid.NY):
                for i in range(CGrid.NX):
                    indices = np.where(P0 == i + j * CGrid.NX + k * CGrid.NX * CGrid.NY)
                    #Harm means in local x slices
                    rx = FGrid.Rx[i];
                    ry = FGrid.Ry[j];
                    rz = FGrid.Rz[k];

                    # y only growing local partition vector Plocy
                    Plocy = np.ones(rx * ry * rz, dtype=int)
                    for kkk in range(0, rz):
                        for iii in range(0, rx):
                            Plocy[iii + kkk * rx * ry::rx] = iii + kkk * rx

                    vol_loc=np.bincount(Plocy,weights=vol[indices])
                    KC3loc=vol_loc/np.bincount(Plocy,weights=Kvol[indices])
                    #Arithm mean of all the local (i,j,k) x slices
                    KC[i,j,k] = np.sum(KC3loc*vol_loc)/np.sum(vol_loc)
        self.Coarse_Mod.UpdateCellData(varname=scalar, array=KC.reshape((-1,1),order='F'))
        m=np.min(KC)
        M=np.max(KC)
        plot_hist(self.Coarse_Mod.GRDECL_Data.SpatialDatas[scalar],\
                  varname=scalar + ": " + "harmonicy: min %3f - max %3f" %(m,M))

    def Upscale_Harmz_mean(self,FineModel, scalar):
        KF=FineModel.GRDECL_Data.SpatialDatas[scalar]
        FGrid = FineModel.GRDECL_Data
        CGrid = self.Coarse_Mod.GRDECL_Data
        vol = FGrid.VOL
        P0 = self.Coarse_Mod_Partition
        # Computing coarse K values
        KC = np.zeros((CGrid.NX, CGrid.NY, CGrid.NZ), dtype=float)
        Kvol=vol/KF
        for k in range(CGrid.NZ):
            for j in range(CGrid.NY):
                for i in range(CGrid.NX):
                    indices = np.where(P0 == i + j * CGrid.NX + k * CGrid.NX * CGrid.NY)

                    #Harm means in local x slices
                    # 3c. Harmonic(in z)-Arithmetic Average
                    rx = FGrid.Rx[i];
                    ry = FGrid.Ry[j];
                    rz = FGrid.Rz[k];

                    # z only growing local partition vector Plocy
                    Plocz = np.ones(rx * ry * rz, dtype=int)
                    for iii in range(0, rx * ry):
                        Plocz[iii::rx * ry] = iii

                    vol_loc=np.bincount(Plocz,weights=vol[indices])
                    KC3loc=vol_loc/np.bincount(Plocz,weights=Kvol[indices])
                    #Arithm mean of all the local (i,j,k) x slices
                    KC[i,j,k] = np.sum(KC3loc*vol_loc)/np.sum(vol_loc)
        self.Coarse_Mod.UpdateCellData(varname=scalar, array=KC.reshape((-1,1),order='F'))
        m=np.min(KC)
        M=np.max(KC)
        plot_hist(self.Coarse_Mod.GRDECL_Data.SpatialDatas[scalar],\
                  varname=scalar + ": " + "harmonicz: min %3f - max %3f" %(m,M))
