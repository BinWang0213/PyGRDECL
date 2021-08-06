from utils import *
import os,sys

class Upscaler:
    def __init__(self):
        """
        Fields and methods to upscale a geomodel
        """
        self.FineModel=None
        self.Coarse_Mod=None
        self.Coarse_Mod_Partition=[]
        self.local_Mod=None
        self.local_Mod_extended=None
        self.nlayer=0
        self.ilayer0=[]
        self.localsize0=[]
        self.Glob_ind=[]


    def create_coarse_model(self):
        fineFname=self.FineModel.fname
        self.Coarse_Mod.fname = os.path.splitext(fineFname)[0] + '_Coarse' +  os.path.splitext(fineFname)[1]
        self.Coarse_Mod_Partition = self.FineModel.GRDECL_Data.fill_coarse_grid(self.Coarse_Mod)
        self.Upscale_Arithmetic_mean("PORO")
        return self.Coarse_Mod

    def create_local_model(self,ind):
        CGrid=self.Coarse_Mod.GRDECL_Data
        assert (ind<CGrid.N),"Ind %Ddlocal exceeds Nx*Ny*Nz: %d"%(ind,CGrid.N)
        fineFname=self.FineModel.fname
        stringloc='_Loc_'+str(ind)+"_nlayer_"+str(self.nlayer)
        self.local_Mod.fname = os.path.splitext(fineFname)[0] +stringloc+  os.path.splitext(fineFname)[1]
        self.Glob_ind,localsize=self.FineModel.GRDECL_Data.get_partition_indices(CGrid,self.nlayer,ind)
        Glob_ind0, self.localsize0 = self.FineModel.GRDECL_Data.get_partition_indices(CGrid, nlayer=0, ind=ind)
        self.ilayer0 = []
        for i, ind in enumerate(self.Glob_ind):
            if ind in Glob_ind0:
                self.ilayer0.append(i)
        self.FineModel.GRDECL_Data.fill_local_grid(self.local_Mod,self.Coarse_Mod_Partition,ind,self.Glob_ind,localsize)
        return self.local_Mod

    def list_upscaling_methods(self):
        mean_methods = ['Arithmetic', 'Geometric', 'Harmonic', 'Harmx', 'Harmy', 'Harmz']
        upsc_methods = [s + "_mean" for s in mean_methods]

        tpfa_methods = ['loc', 'glob', 'loc_vol_average']
        upsc_methods += ["TPFA_" + s for s in tpfa_methods]
        return upsc_methods

    def Upscale_Arithmetic_mean(self, scalar):
        FGrid = self.FineModel.GRDECL_Data
        Finefield = FGrid.SpatialDatas[scalar]
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

    def Upscale_Geometric_mean(self, scalar):
        FGrid = self.FineModel.GRDECL_Data
        Finefield = FGrid.SpatialDatas[scalar]
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

    def Upscale_Harmonic_mean(self, scalar):
        FGrid = self.FineModel.GRDECL_Data
        Finefield = FGrid.SpatialDatas[scalar]
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

    def Upscale_Harmx_mean(self, scalar):
        FGrid = self.FineModel.GRDECL_Data
        KF = FGrid.SpatialDatas[scalar]
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

    def Upscale_Harmy_mean(self, scalar):
        FGrid = self.FineModel.GRDECL_Data
        KF = FGrid.SpatialDatas[scalar]
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

    def Upscale_Harmz_mean(self, scalar):
        FGrid = self.FineModel.GRDECL_Data
        KF = FGrid.SpatialDatas[scalar]
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

    def Upscale_TPFA_loc(self):
        CGrid = self.Coarse_Mod.GRDECL_Data
        print("[UPSCALING TPFA loc], nlayer:",self.nlayer)
        for ind in range(CGrid.N):
            sys.stdout.write("Local computations for coarse cell: %d / %d \r" % (ind+1,CGrid.N))
            sys.stdout.flush()
            # print("local upscaling for block: ",ind)
            Model3= self.FineModel.create_local_model(ind)
            Nx,Ny,Nz=self.localsize0
            N  = Nx*Ny*Nz
            dx = np.array(Model3.GRDECL_Data.DX)[self.ilayer0]
            dy = np.array(Model3.GRDECL_Data.DY)[self.ilayer0]
            dz = np.array(Model3.GRDECL_Data.DZ)[self.ilayer0]

            directions = ["i", "j", "k"]
            Lx = np.sum(dx[:Nx])
            Ly = np.sum(dy[0:Nx * Ny:Nx])
            Lz = np.sum(dz[0:N:Nx * Ny])
            Lxyz = [Lx, Ly, Lz]
            dxyz = [dx, dy, dz]
            scalars = ["PERMX", "PERMY", "PERMZ"]
            for i, dir in enumerate(directions):
                Model3.compute_TPFA_Pressure(Press_inj=1, direction=dir)
                GradP, V = Model3.computeGradP_V()
                GradP =np.array( GradP.reshape((3, Model3.GRDECL_Data.N), order='F'))[:,self.ilayer0]
                V = np.array(V.reshape((3, Model3.GRDECL_Data.N), order='F'))[:,self.ilayer0]

                Ind_face0, Ind_face1 = Model3.compute_bdry_indices(direction=dir)
                Ind_face0=[ind for ind in Ind_face0 if ind in self.Glob_ind[self.ilayer0] ]
                Ind_face1=[ind for ind in Ind_face1 if ind in self.Glob_ind[self.ilayer0] ]

                L = Lxyz[i]
                A = Lxyz[(i + 1) % 3] * Lxyz[(i + 2) % 3]
                q = np.sum(
                    V[i, Ind_face1] * np.array(dxyz[(i + 1) % 3])[Ind_face1] * np.array(dxyz[(i + 2) % 3])[Ind_face1])
                CGrid.SpatialDatas[scalars[i]][ind] = np.abs(q) * L / A
        i = 0
        m = np.min(CGrid.SpatialDatas[scalars[i]])
        M = np.max(CGrid.SpatialDatas[scalars[i]])
        plot_hist(CGrid.SpatialDatas[scalars[i]], \
                  varname=scalars[i] + ": " + "Local Flow min %3f - max %3f" % (m, M))

    def Upscale_TPFA_glob(self):
        FGrid=self.FineModel.GRDECL_Data
        CGrid = self.Coarse_Mod.GRDECL_Data
        scalars = ["PERMX", "PERMY", "PERMZ"]
        directions = ["i", "j", "k"]
        for i, dir in enumerate(directions):
            self.FineModel.compute_TPFA_Pressure(Press_inj=1, direction=dir)
            GradP, V = self.FineModel.computeGradP_V()
            V = V.reshape((3, FGrid.N),order='F')
            GradP = GradP.reshape((3, FGrid.N),order='F')

            for ind_model_local in range(CGrid.N):
                Model3=self.FineModel.create_local_model(ind_model_local)
                # Model3 = self.local_Mod
                V_loc=np.array(V,dtype=float)[:,self.Glob_ind]
                GradP_loc=np.array(GradP,dtype=float)[0:3,self.Glob_ind]

                Nx = Model3.GRDECL_Data.NX
                Ny = Model3.GRDECL_Data.NY
                Nz = Model3.GRDECL_Data.NZ
                N = Model3.GRDECL_Data.N
                dx = Model3.GRDECL_Data.DX
                dy = Model3.GRDECL_Data.DY
                dz = Model3.GRDECL_Data.DZ

                dxyz = [dx, dy, dz]
                Ind_face0, Ind_face1 = Model3.compute_bdry_indices(direction=dir)
                GradP_face = np.sum(GradP_loc[i, Ind_face1] * \
                           np.array(dxyz[(i + 1) % 3])[Ind_face1] *\
                           np.array(dxyz[(i + 2) % 3])[Ind_face1])
                q = np.sum(V_loc[i, Ind_face1] * \
                           np.array(dxyz[(i + 1) % 3])[Ind_face1] *\
                           np.array(dxyz[(i + 2) % 3])[Ind_face1])
                CGrid.SpatialDatas[scalars[i]][ind_model_local] = np.abs(q)  / (GradP_face)
        i=0
        m=np.min(self.Coarse_Mod.GRDECL_Data.SpatialDatas[scalars[i]])
        M=np.max(self.Coarse_Mod.GRDECL_Data.SpatialDatas[scalars[i]])
        plot_hist(self.Coarse_Mod.GRDECL_Data.SpatialDatas[scalars[i]], \
                  varname=scalars[i] + ": " + "Global Flow min %3f - max %3f" %(m,M))

    def Upscale_TPFA_loc_vol_average(self):
        FGrid=self.FineModel.GRDECL_Data
        CGrid = self.Coarse_Mod.GRDECL_Data
        print("[UPSCALING TPFA loc volume average], nlayer:",self.nlayer)

        for ind in range(CGrid.N):
            sys.stdout.write("Local computations for coarse cell: %d / %d \r" % (ind+1,CGrid.N))
            sys.stdout.flush()
            Model3=self.FineModel.create_local_model(ind)
            # Model3 = self.local_Mod
            Nx,Ny,Nz=self.localsize0
            N  = Nx*Ny*Nz
            dx = np.array(Model3.GRDECL_Data.DX)[self.ilayer0]
            dy = np.array(Model3.GRDECL_Data.DY)[self.ilayer0]
            dz = np.array(Model3.GRDECL_Data.DZ)[self.ilayer0]

            scalars = ["PERMX","PERMXY", "PERMXZ","PERMYX", "PERMY","PERMYZ","PERMZX","PERMZY","PERMZZ"]
            directions = ["i", "j", "k"]
            # Initialize local system to be solved
            A=np.zeros((9+3,9))
            # Enforce symmetry
            A[9,1]=1;A[9,3]=-1
            A[10,2]=1;A[10,6]=-1
            A[11,5]=1;A[11,7]=-1

            b=np.zeros((9+3,1))
            vol=dx*dy*dz
            vol_tot = np.sum(vol)

            for ind_sys,dir in enumerate(directions):
                Model3.compute_TPFA_Pressure(Press_inj=1, direction=dir)
                GradP, V = Model3.computeGradP_V()
                GradP =np.array( GradP.reshape((3, Model3.GRDECL_Data.N), order='F'))[:,self.ilayer0]
                V = np.array(V.reshape((3, Model3.GRDECL_Data.N), order='F'))[:,self.ilayer0]
                for ivar in range(3):
                    Mean_GradPi = np.sum(GradP[ivar, :] *vol) /vol_tot
                    Mean_Vi = np.sum(V[ivar, :] *vol) /vol_tot
                    irow=ivar+3*ind_sys
                    b[irow]       = -Mean_Vi
                    A[0+3*ind_sys,ivar]  = Mean_GradPi
                    A[1+3*ind_sys,ivar+3]= Mean_GradPi
                    A[2+3*ind_sys,ivar+6]= Mean_GradPi
            K_loc, residuals, rank, s=np.linalg.lstsq(A,b)
            for iscalar,scalar in enumerate(scalars):
                if scalar not in CGrid.SpatialDatas:
                    self.Coarse_Mod.CreateCellData(varname=scalar, val=1)
                CGrid.SpatialDatas[scalar][ind] = K_loc[iscalar]
        i=0
        m=np.min(CGrid.SpatialDatas[scalars[i]])
        M=np.max(CGrid.SpatialDatas[scalars[i]])
        plot_hist(CGrid.SpatialDatas[scalars[i]], \
                  varname=scalars[i] + ": " + "Local Flow Vol Avg min %3f - max %3f" %(m,M))

    def Upscale_Perm(self, upsc_methods):
        if isinstance(upsc_methods, str):
            upsc_methods = [upsc_methods]

        scalars = ["PERMX", "PERMY", "PERMZ"]
        if upsc_methods[0][-4:] == "mean":
            if (len(upsc_methods)==1): upsc_methods *=3
            for i, upsc_method in enumerate(upsc_methods):
                method = "Upscale_" + upsc_method
                scalar = scalars[i]
                getattr(self, method)(scalar)
        else:
            getattr(self, "Upscale_" + upsc_methods[0])()