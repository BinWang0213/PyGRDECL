#########################################################################
#       (C) 2017-2018 Department of Petroleum Engineering,              # 
#       Univeristy of Louisiana at Lafayette, Lafayette, US.            #
#                                                                       #
# This code is released under the terms of the BSD license, and thus    #
# free for commercial and research use. Feel free to use the code into  #
# your own project with a PROPER REFERENCE.                             #
#                                                                       #
# PyGRDECL Code                                                         #
# Author: Bin Wang                                                      # 
# Email: binwang.0213@gmail.com                                         # 
#########################################################################

import os
import numpy as np
import math
SupportKeyWords=[
    'SPECGRID', #Dimenion of the corner point grid
    'DIMENS',   #Define the dimension of the cartesian grid
    'TOPS','DX','DY','DZ',
    'COORD','ZCORN','ACTNUM'
    'PORO',
    'PERMX' , 'PERMXY', 'PERMXZ', 
    'PERMYX', 'PERMY' , 'PERMYZ', 
    'PERMZX', 'PERMZY', 'PERMZ',
    'SW_NPSL'
]

KeyWordsDatatypes=[#Corrsponding data types
    int,
    float,
    int,int,int,int,
    float,float,
    float,
    float,float,float,
    float,float,float,
    float,float,float
]

class nodes:
    def __init__(self):
        self.coords=[]
        self.num=0

class cells:
    def __init__(self):
        self.num=0
        self.indexMap=[]

class faces:
    def __init__(self):
        self.nodes    = np.empty((0,1), int)
        self.nodePos  = np.empty((0,1), int)
        self.neighbors= np.empty((0,1), int)
        self.tag=[]
        self.cellTags=[]

class GRDECL_Parser:
    def __init__(self,filename='',nx=0,ny=0,nz=0):
        """Eclipse Input file(GRDECL) Parser 
        Keywords Reference: file format:http://petrofaq.org/wiki/Eclipse_Input_Data

        Arguments
        ---------
        NX, NY, NZ         -- Grid dimension.
        Trans(i01,j01,k01) -- Transmissibility in i,j,k direction
        fault(i01,j01)     -- Fault indicator in i,j direction(0-sealing, 0.5-partially connecting, 1-fully connecting)
        GRID_type           - 0-Cartesian 1-Corner point

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2017
        """
        self.fname=filename
        self.NX=nx
        self.NY=ny
        self.NZ=nz
        self.N=nx*ny*nz
        self.GRID_type='NaN'

        #Cartesian gridblock data KeyWords
        self.TOPS=[]
        self.DX=[]
        self.DY=[]
        self.DZ=[]

        #Corner point gridblock data KeyWrods (not support now)
        self.COORD=[]
        self.ZCORN=[] #http://maoxp9.blog.163.com/blog/static/122653420093894133671/
        self.ACTNUM=[]

        #MRST Format grid attributes
        self.gridDim=3
        self.cartDims=[]
        self.nodes=nodes()
        self.cells=cells()
        self.faces=faces()

        #Petrophysics data Keywords
        self.SpatialDatas={}
        self.SkipedKeywords=0

        #Read GRDECL file when initializing the class
        if(len(filename)>0):
            self.read_GRDECL()

        #Derived variabls 
        self.CELL_FAULT=[]
    
    ######[read_GRDECL]######
    def read_GRDECL(self):
        """Read input file(GRDECL) of Reservoir Simulator- Petrel (Eclipse)  
        file format:http://petrofaq.org/wiki/Eclipse_Input_Data
        
        Arguments
        ---------
        NX, NY, NZ -- Grid dimension.
        blockData_raw -- [0] Keywords [1] values
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2017
        """
        debug=0

        print('[Input] Reading ECLIPSE/PETREL file \"%s\" ....'%(self.fname))

        #Read whole file into list
        f=open(self.fname)
        contents=f.read()
        contents=RemoveCommentLines(contents,commenter='--')
        contents_in_block=contents.strip().split('/') #Sepeart input file by slash /
        contents_in_block = [x for x in contents_in_block if x]#Remove empty block at the end
        NumKeywords=len(contents_in_block)

        GoodFlag=0
        for i,block in enumerate(contents_in_block):#Keyword, Block-wise
            #Clean the data where no spliter \ provided
            block=scanKeyword(block)

            blockData_raw=block.strip().split()
            Keyword=''
            DataArray=[]
            if(len(blockData_raw)>1):
                if(blockData_raw[0]=='ECHO'): #This keyword may next to real keyword
                    Keyword,DataArray=blockData_raw[1],blockData_raw[2:]                    
                else:
                    Keyword,DataArray=blockData_raw[0],blockData_raw[1:]

            #Read Grid Dimension [SPECGRID] or [DIMENS] 
            if(Keyword=='DIMENS'):
                DataArray=np.array(DataArray[:3],dtype=int)
                self.GRID_type='Cartesian'
                self.NX,self.NY,self.NZ=DataArray[0],DataArray[1],DataArray[2]
                self.N=self.NX*self.NY*self.NZ
                print("     Grid Type=%s Grid" %(self.GRID_type))
                print("     Grid Dimension(NX,NY,NZ): (%s x %s x %s)"%(self.NX,self.NY,self.NZ))
                print("     NumOfGrids=%s"%(self.N))
                print('     NumOfKeywords=%s'%(NumKeywords))
                print("     Reading Keyword %d [%s] " %(i+1,Keyword),end='')
                GoodFlag=1
                continue
            elif(Keyword=='SPECGRID'):
                DataArray=np.array(DataArray[:3],dtype=int)
                self.GRID_type='CornerPoint'
                self.NX,self.NY,self.NZ=DataArray[0],DataArray[1],DataArray[2]
                self.N=self.NX*self.NY*self.NZ
                print("     Grid Type=%s" %(self.GRID_type))
                print("     Grid Dimension(NX,NY,NZ): (%s x %s x %s)"%(self.NX,self.NY,self.NZ))
                print("     NumOfGrids=%s"%(self.N))
                print('     NumOfKeywords=%s'%(NumKeywords))
                print("     Reading Keywords [%s] " %(Keyword),end='')
                GoodFlag=1
                continue
            
            if(self.GRID_type=='NaN'):#Skip unnecessary keywords
                continue

            if(Keyword in SupportKeyWords): #We need parse the special format in 
                if(len(DataArray)==1 and '.' in DataArray[0]):
                    folder_name=os.path.dirname(self.fname)
                    DataArray=self.read_IncludeFile(os.path.join(folder_name,DataArray[0]),self.N)
                #print(Keyword,DataArray)
                DataArray=parseDataArray(DataArray)
            

            #Read Grid spatial information, x,y,z ordering
            #Corner point cell
            if(Keyword=='COORD'):# Pillar coords
                assert len(DataArray)==6*(self.NX+1)*(self.NY+1),'[Error] Incompatible COORD data size!'
                self.COORD=np.array(DataArray,dtype=float)       
            elif(Keyword=='ZCORN'):# Depth coords
                assert len(DataArray)==8*self.N, '[Error] Incompatible ZCORN data size!'
                self.ZCORN=np.array(DataArray,dtype=float)
            
            #Cartesian cell
            elif(Keyword=='DX'):# Grid size in X dir
                assert len(DataArray)==self.N, '[Error] Incompatible DX data size!'
                self.DX=np.array(DataArray,dtype=float)
            elif(Keyword=='DY'):# Grid size in Y dir
                assert len(DataArray)==self.N, '[Error] Incompatible DY data size!'
                self.DY=np.array(DataArray,dtype=float)
            elif(Keyword=='DZ'):# Grid size in Z dir
                assert len(DataArray)==self.N, '[Error] Incompatible DZ data size!'
                self.DZ=np.array(DataArray,dtype=float)
            elif(Keyword=='TOPS'):# TOP position
                assert len(DataArray)==self.N, '[Error] Incompatible TOPS data size!'
                self.TOPS=np.array(DataArray,dtype=float)

            #Read Grid Properties information
            else:
                self.LoadVar(Keyword,DataArray,DataSize=self.N)

        f.close()
        assert GoodFlag==1,'Can not find grid dimension info, [SPECGRID] or [DIMENS]!'
        print('.....Done!')


        #Genetrate TOPS for cartesian grid if TOPS if not given
        if(self.GRID_type=='Cartesian' and len(self.TOPS)==0):
            self.TOPS=np.zeros(self.N)
            for k in range(self.NZ-1):
                for j in range(self.NY):
                    for i in range(self.NX):
                        ijk=getIJK(i,j,k,self.NX,self.NY,self.NZ)
                        ijk_next=getIJK(i,j,k+1,self.NX,self.NY,self.NZ)
                        self.TOPS[ijk_next] = self.TOPS[ijk] + self.DZ[ijk]
    

    def buildCartGrid(self,physDim=[100.0,100.0,10.0],gridDims=[10,10,1]):
        """Build simple cartesian grid 
        
        Arguments
        ---------
        physDim  -- physical dimensions of system
        gridDims -- grid dimension of system

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Feb. 2019
        """
        self.NX,self.NY,self.NZ=gridDims
        self.N=self.NX*self.NY*self.NZ
        self.GRID_type='Cartesian'


        print("     Grid Type=%s Grid" %(self.GRID_type))
        print("     Grid Dimension(NX,NY,NZ): (%s x %s x %s)"%(self.NX,self.NY,self.NZ))
        print("     NumOfGrids=%s"%(self.N))

        #Assign value to cart grid
        self.DX=np.ones(self.N)*physDim[0]/self.NX
        self.DY=np.ones(self.N)*physDim[1]/self.NY
        self.DZ=np.ones(self.N)*physDim[2]/self.NZ
        self.TOPS=np.zeros(self.N)
        for k in range(self.NZ-1):
                for j in range(self.NY):
                    for i in range(self.NX):
                        ijk=getIJK(i,j,k,self.NX,self.NY,self.NZ)
                        ijk_next=getIJK(i,j,k+1,self.NX,self.NY,self.NZ)
                        self.TOPS[ijk_next] = self.TOPS[ijk] + self.DZ[ijk]
        
        #Build up basic spatial propertis
        self.SpatialDatas["PERMX"]=np.ones(self.N)*10.0
        self.SpatialDatas["PERMY"]=np.ones(self.N)*10.0
        self.SpatialDatas["PERMZ"]=np.ones(self.N)*10.0
        self.SpatialDatas["PORO"]=np.ones(self.N)*0.3


    def buildCPGGrid(self, physDims=[1.0, 1.0, 0.5], gridDims=[3, 3, 3],opt={'disturbed':True,'flat':False},faultDrop=0.):
        """Build simple corner point grid
        Arguments
        ---------
        physDim   -- physical dimensions of system
        gridDims  -- grid dimension of system
        opt['disturbed'] -- disturb from cartesian grid: skew pillars
        opt['flat'] -- sinusoidal z
        faultDrop -- fault drop at half the CPG in x direction


        Author:Mustapha Zakari(mustapha.zakari@univ-lorraine.fr)
        Date: May. 2021
        Python version based on a translation of MRST simpleGrdecl
        """
        self.NX, self.NY, self.NZ = gridDims
        self.N = self.NX * self.NY * self.NZ
        self.GRID_type = 'CornerPoint'

        print("     Grid Type=%s Grid" % (self.GRID_type))
        print("     Grid Dimension(NX,NY,NZ): (%s x %s x %s)" % (self.NX, self.NY, self.NZ))
        print("     NumOfGrids=%s" % (self.N))

        # Create x,y,z vectors. Disturb x,z
        x_= np.linspace(0, physDims[0], self.NX + 1)
        y_= np.linspace(0, physDims[1], self.NY + 1)
        z_= np.linspace(0, physDims[2], self.NZ + 1)
        x, y, z = np.meshgrid(x_, y_, z_, indexing='ij')

        if (opt['disturbed']):
            # skew pillars
            xmid=x[self.NX//2,0,0];
            zmid=z[0,0,self.NZ//2];
            x = x +3*0.05* (xmid - abs(x - xmid)) * (z - physDims[2])/physDims[2];#disturb more nead xmid and zmax
            x = physDims[0] * x/x.max()#rescale
        if (not opt['flat']):
            # add sinusoidal variations for z
            def fz(x, y):
                z =  np.sin(math.pi * x) + 0.5*np.sin(math.pi * (y +2*x));
                return z
            for k in range(gridDims[2]+1):
                xi = x[:, :, k] / physDims[0];
                eta = y[:, :, k] / physDims[1];
                z[:, :, k] = z[:, :, k] - 0.05*physDims[0]*fz(xi, eta);
            z = np.sort(z, 2);
            z = physDims[2]* (z-z.min())/(z.max()-z.min())


        # Create pillars
        npillar = (self.NX + 1)*(self.NY + 1)
        ## Each pillar bottom & top xyz COORD (6 coordinates)
        lines = np.zeros((npillar, 6));
        ## Set columns Bottom/Top (z[0]/z[end]) coordinates: xB yB zB xT yT zT
        for i, V in enumerate([x, y, z]):
            lines[:, i + 0] = V[:, :,  0].reshape((npillar), order='F');      #Bottom:iz=0
            lines[:, i + 3] = V[:, :, -1].reshape((npillar), order='F'); #Top:iz=lastindex=-1
        ## flatten lines to create COORD
        self.COORD = lines.reshape((np.prod(lines.size)))

        # Assign z-coordinates
        ## Add repeated indices to z in each direction:
        ## 0,    1,1,       2,2,   ... ,gridDims(d)-1,gridDims(d)-1, gridDims(d)
        ## 0, 2//2,3//2, 4//2,5//2,... ,                           ,2*gridDims(d)//2
        def repInd(d):
            return np.linspace(1, 2*gridDims[d], 2*gridDims[d],dtype=int)//2
        IX, IY, IZ = repInd(0), repInd(1), repInd(2)

        self.ZCORN = np.zeros((2 * self.NX, 2 * self.NY, 2 * self.NZ), dtype=float)
        for k in range(2 * self.NZ):
            for j in range(2 * self.NY):
                for i in range(2 * self.NX):
                    self.ZCORN[i, j, k] = z[IX[i], IY[j], IZ[k]]

        # Add a fault drop after half the x range: (odd index)
        ## [0, 1],[2, 3],4, ... ,2*(gridDims(d)-1)],[2*gridDims(d)-1, 2*gridDims(d)]
        iFault=self.NX//2
        self.ZCORN[2*iFault:,:,:]+=faultDrop

        # flatten ZCORN
        self.ZCORN = self.ZCORN.reshape((np.prod(self.ZCORN.size)), order='F')
        # Build up basic spatial propertis
        self.ACTNUM = np.ones((self.N))
        self.SpatialDatas["PERMX"] = np.ones(self.N) * 10.0
        self.SpatialDatas["PERMY"] = np.ones(self.N) * 10.0
        self.SpatialDatas["PERMZ"] = np.ones(self.N) * 10.0
        self.SpatialDatas["PORO"] = np.ones(self.N) * 0.3

    def buildCornerPointNodes(self):
        # Construct nodal coordinates for Corner Point Grid
        nx, ny, nz = self.NX, self.NY, self.NZ
        ncell=nx*ny*nz
        npillar = (nx + 1) * (ny + 1)

        # Enumerate pillars in grid associate one pillar per node
        # Id from  reshaping of [1 2 3 ... npillar] (column major Fortran ordering 11,21,12,22)
        pillarId = (np.linspace(1, npillar, npillar,dtype=int)).reshape((nx + 1, ny + 1), order='F')
        pillarId-=1

        # Assign separately 4 pillar Ids per cell 00,10,01,11
        p1 = pillarId[0:nx  ,0:ny];
        p2 = pillarId[1:nx+1,0:ny];
        p3 = pillarId[0:nx  ,1:ny+1];
        p4 = pillarId[1:nx+1,1:ny+1];
        del pillarId;
        # Assign 4 pillar Ids per cell (with repetitions) and vector reshape it
        lineID=np.zeros((2*nx,2*ny,2*nz),dtype=int)
        for k in range(2*nz):
            lineID[0:2*nx:2, 0:2*ny:2,k] = p1;
            lineID[1:2*nx:2, 0:2*ny:2,k] = p2;
            lineID[0:2*nx:2, 1:2*ny:2,k] = p3;
            lineID[1:2*nx:2, 1:2*ny:2,k] = p4;
            lineID2=lineID.reshape((8*ncell),order='F')
        del lineID,p1,p2,p3,p4;

        # Use COORD and 4 pillar Ids per cell to get
        # xmymzmxMyMzM pillar coord for their 8 cell nodes
        COORD2=self.COORD.reshape((npillar,6))
        node_pillar=np.zeros((8*ncell,6),dtype=float)
        for k in range(8*ncell):
            for j in range(6):
                node_pillar[k,j]=COORD2[lineID2[k],j]
        del COORD2,lineID2;

        # Reconstruct nodal coordinates from pillars and ZCORN
        linFactor = (self.ZCORN[:] - node_pillar[:, 2]) / (node_pillar[:, 5] - node_pillar[:, 2]);
        self.X=node_pillar[:,0]+linFactor*(node_pillar[:,3]-node_pillar[:,0])
        self.Y=node_pillar[:,1]+linFactor*(node_pillar[:,4]-node_pillar[:,1])

        self.X=(self.X).reshape((2*nx,2*ny,2*nz),order='F')
        self.Y=(self.Y).reshape((2*nx,2*ny,2*nz),order='F')
        self.Z=(self.ZCORN).reshape((2*nx,2*ny,2*nz),order='F')

        ## Reverse z if dz<0
        # Expand actnum by 2 in each direction
        # from numpy import tile
        # self.ACTNUM = (self.ACTNUM).reshape((self.NX, self.NY, self.NZ), order='F')
        # a = tile(self.ACTNUM, (2, 2, 2))
        # z=self.Z; z[a==0]=float('NaN')
        # dz=np.diff(z)
        # self.Z=dz;

        # Add top+bottom layers to ensure correct processing of outer bdry at faults
        minz = (self.Z).min();
        maxz = (self.Z).max();

        e = np.zeros((2 * self.NX, 2 * self.NY, 1))
        self.Z = np.concatenate((minz - 2 + e, minz - 1 + e, self.Z, maxz + 1 + e, maxz + 2 + e), axis=2)

        e1 = (self.X[:, :, 0]).reshape((2 * self.NX, 2 * self.NY, 1));
        e2 = (self.X[:, :, -1]).reshape((2 * self.NX, 2 * self.NY, 1));
        self.X = np.concatenate((e1, e1, self.X, e2, e2), axis=2)

        e1 = (self.Y[:, :, 0]).reshape((2 * self.NX, 2 * self.NY, 1));
        e2 = (self.Y[:, :, -1]).reshape((2 * self.NX, 2 * self.NY, 1));
        self.Y = np.concatenate((e1, e1, self.Y, e2, e2), axis=2)

        # Mark active new layers
        actnum = (self.ACTNUM).reshape(self.NX, self.NY, self.NZ);
        e = np.ones((self.NX, self.NY, 1), dtype=bool)
        self.actnum = np.concatenate((e, actnum, e), axis=2)
        self.numAuxiliaryCells = 2 * np.prod(e.shape)

        # Replace nan coordinates in X,Y,Z by inf to avoid nan
        from numpy import isnan
        for V in [self.X, self.Y, self.Z]:
            V[isnan(V)] = float('Inf')

    def createInitialGrid(self):
        # Find unique points
        self.GRID_type='INVALID'
        self.dim=3;
        N=np.prod(self.Z.size);
        # Return the sorted list/array of unique elements of the coord array
        # sorted following i,j,k implies Z must be the first lexicographic axis
        self.nodes.coords=np.zeros((N,3))
        self.nodes.coords[:,0]=self.Z.reshape((N), order='F')
        self.nodes.coords[:,1]=self.Y.reshape((N), order='F')
        self.nodes.coords[:,2]=self.X.reshape((N), order='F')

        # ia: indices of kept unique values in initial array: a,ia=unique(ar0) <=> a=ar0[ia]
        # ic: indices to reconstruct initial array: a,ia,ic=unique(ar0) <=> ar0=a[ic]
        self.nodes.coords, ia,ic =np.unique(self.nodes.coords, return_index=True, return_inverse=True, axis=0)
        self.nodes.coords=np.fliplr( self.nodes.coords)
        self.nodes.num=self.nodes.coords[:,0].size
        self.cartDims= [number // 2 for number in self.X.shape] #nb of cells in each direction

        self.cells.num=np.prod(self.cartDims)
        # Cell indices growing i,j,k
        self.cells.indexMap=np.array(range(self.cells.num))

        # Store indices in P to reconstruct and repeat all point coords
        self.P=ic.reshape(self.X.shape,order='F')
        # Assign each cart cell each unique Id in array B : growing i,j,k
        self.B=np.array(range(self.cells.num)).reshape(self.cartDims,order='F')

    def index(self,i,j,k,sz):
        I,J,K=  np.meshgrid(i,j,k, indexing='ij')
        #find linear index from 3D indices equiv to sub2ind in matlab
        return np.ravel_multi_index([I,J,K], sz, order='F').reshape(-1,order='F')

    def findFaces(self):
        # Find regular faces including their nodes, cell neighbors and direction
        #         PARAMETERS:
        #    G      - Grid struct that faces will be filled into.
        #
        #    P      - 3d-array of point numbers for each cell in the grid, .i.e.,
        #             the grid cell (i,j,k) has corner points numbers
        #                 n1 = P(2*i-1,2*j-1,2*k-1)
        #                 n2 = P(2*i  ,2*j-1,2*k-1)
        #                 n3 = P(2*i-1,2*j  ,2*k-1)
        #                 n4 = P(2*i  ,2*j  ,2*k-1)
        #                 n5 = P(2*i-1,2*j-1,2*k  )
        #                 n6 = P(2*i  ,2*j-1,2*k  )
        #                 n6 = P(2*i-1,2*j  ,2*k  )
        #                 n8 = P(2*i  ,2*j  ,2*k  )
        #
        #    B      - Block numbers in the grid.  These are included for convenience
        #             as they could easily have been computed from (i,j,k).
        #
        #    actnum - Array of 0/1, 0 indicate inactive cell, 1 active cell.
        #
        #    tags   - Face tag (numeric) that should be inserted into G.
        #
        #    opt    - options struct from which opt.verbose is used.
        #
        sz  = self.P.shape
        szb = self.B.shape

        # Find face corners
        di=1;#i step next index
        dj=sz[0]#j step next index
        dk=dj*sz[1]#k step next index
        # first internal faces (West)
        k=self.index(np.array(range(1,sz[0]-1,2)),np.array(range(0,sz[1],2)),np.array(range(0,sz[2]-1,2)),sz)
        # Internal faces
        # converting unique indices k in first global indices from P
        # f = [P(k), P(k+dj), P(k+dj+dk), P(k+dk)];
        #     [c1    c2       c3          c4     ]
        c1=self.P[np.unravel_index(k, self.P.shape, 'F')]
        c2=self.P[np.unravel_index(k+dj, self.P.shape, 'F')]
        c3=self.P[np.unravel_index(k+dj+dk, self.P.shape, 'F')]
        c4=self.P[np.unravel_index(k+dk, self.P.shape, 'F')]
        self.f=np.ones((c1.size,4))
        C=[c1,c2,c3,c4]
        for i in range(4):
            self.f[:,i]=C[i]
        # first internal unique indices of faces (East)
        k=k+di
        self.g=np.ones((c1.size,4))
        c1=self.P[np.unravel_index(k, self.P.shape, 'F')]
        c2=self.P[np.unravel_index(k+dj, self.P.shape, 'F')]
        c3=self.P[np.unravel_index(k+dj+dk, self.P.shape, 'F')]
        c4=self.P[np.unravel_index(k+dk, self.P.shape, 'F')]
        C=[c1,c2,c3,c4]
        for i in range(4):
            self.g[:,i]=C[i]

        # Neighbor cells
        k=self.index(np.array(range(szb[0]-1)),np.array(range(szb[1])),np.array(range(szb[2])),szb)
        self.c1 = self.B[np.unravel_index(k, self.B.shape, 'F')]
        k=self.index(np.array(range(1,szb[0])),np.array(range(szb[1])),np.array(range(szb[2])),szb)
        self.c2 = self.B[np.unravel_index(k, self.B.shape, 'F')]

        # Keep only perfectly matching cell pairs (c1,c2)
        self.h=np.all(self.f==self.g,1)
        shape=[i//2 for i in sz];shape[0]-=1
        self.h=np.all(self.h.reshape((shape),order='F'),2)
        self.h = np.repeat(self.h[:, :, np.newaxis], sz[2]//2, axis=2)
        self.c1=self.c1[np.where(self.h.reshape((self.h.size),order='F'))]
        self.c2=self.c2[np.where(self.h.reshape((self.h.size),order='F'))]
        self.h=np.where(self.h.reshape((self.h.size), order='F'))
        self.f = self.f[self.h[0][:],:]
        del self.g,self.h

        # Regular Boundary faces
        k=self.index(np.array([0,sz[0]-1]),np.array(range(0,sz[1],2)),np.array(range(0,sz[2]-1,2)),sz)
        # fB = [self.P(k), P(k + dj), P(k + dj + dk), P(k + dk)];
        c1=self.P[np.unravel_index(k, self.P.shape, 'F')]
        c2=self.P[np.unravel_index(k+dj, self.P.shape, 'F')]
        c3=self.P[np.unravel_index(k+dj+dk, self.P.shape, 'F')]
        c4=self.P[np.unravel_index(k+dk, self.P.shape, 'F')]
        C=[c1,c2,c3,c4]
        self.fB=np.ones((c1.size,4))
        for i in range(4):
            self.fB[:,i]=C[i]
        k=self.index([0,szb[0]-1],np.array(range(szb[1])),np.array(range(szb[2])),szb)
        self.cB = self.B[np.unravel_index(k, self.B.shape, 'F')]
        # # tagB = repmat(tags(:), [numel(cB) / 2, 1]);
        self.tagB= np.tile(self.tag, self.cB.size//2)

        # Append boundary faces
        self.f=np.vstack((self.f,self.fB))
        del self.fB
        self.cB[self.tagB==self.tag[0]]=-1
        self.c1=np.hstack((self.c1,self.cB))#
        self.cB = self.B[np.unravel_index(k, self.B.shape, 'F')]
        self.cB[self.tagB==self.tag[1]]=-1
        self.c2=np.hstack((self.c2,self.cB))#

        # Filter out inactive and degenerate cell pairs
        # Remove inactive cells
        self.actnum_lin=np.reshape(self.actnum,-1)
        self.c1[self.c1 != -1] *= (self.actnum_lin[self.c1 != -1]).astype('int32')
        self.c2[self.c2 != -1] *= (self.actnum_lin[self.c2 != -1]).astype('int32')

        # Remove faces with no neighbors and pinched faces
        ind = ((self.c1 != -1) | (self.c2 != -1)) & ~((self.f[:, 0] == self.f[:, 3]) & (self.f[:, 1] == self.f[:, 2]))
        self.f = self.f[ind, :]
        self.c1 = self.c1[ind]
        self.c2 = self.c2[ind]

        # Remove zero-area faces
        ind =(self.f[:,0]==self.f[:,3]) & (self.f[:,1]==self.f[:,2]);
        self.f = self.f[~ind, :]
        self.c1 = self.c1[~ind]
        self.c2 = self.c2[~ind]


        # Remove repeated nodes from pinch
        self.f[self.f[:,0]==self.f[:,3], 0]=float("nan");
        self.f[self.f[:,1]==self.f[:,2], 0]=float("nan");
        self.F=np.reshape(np.hstack((self.f,float("inf")*np.ones((self.f.shape[0],1),dtype=float))),(-1,1))
        self.nF=np.where(self.F==float("inf"))
        self.nF=np.diff(np.hstack(([-1],self.nF[0])))-1
        self.F=self.F[self.F!=float("inf")]

        # Write results to grid srtucture
        n=self.c1.size
        print('Found ',n,' new regular faces')
        # self.faces.nodes =np.empty((0,1),dtype=int)
        self.faces.nodes     = [self.faces.nodes,self.F][1]
        self.faces.neighbors = [self.faces.neighbors,np.transpose(np.vstack((self.c1,self.c2)))][1]
        # self.arr1=np.diff(self.faces.nodePos.astype(float),axis=0),
        self.arr=np.hstack((0,self.nF))
        self.faces.nodePos   = np.cumsum(self.arr);
        self.faces.tag        = [self.faces.tag,        np.zeros((n, 1))][1];
        self.faces.cellTags   = [self.faces.cellTags,   np.tile(self.tag, (self.c1.size, 1))][1];


    def findFaults(self):

        print('findFaults')
    def process_pillar_faces(self):
        print('Processing regular faces')
        self.findFaces()


    def processGRDECL(self):
        # Construct nodal coordinates  from pillars and ZCORN
        # -> self.X,self.Y,self.Z
        self.buildCornerPointNodes()
        # Init Grid nodes,cells attributes
        # Store indices in P to reconstruct all repeated point coords from unique points
        # Store each cart Block/Cell unique Id in array B : lexico growing i,j,k
        # -> self.nodes,self.cells,self.P,self.B
        self.createInitialGrid()

        # Free X,Y,Z spaces
        self.X,self.Y,self.Z=[],[],[]
 
        # Process faces with constant i-index
        self.tag=[1,2] #West East
        self.process_pillar_faces()




    def LoadVar(self,Keyword,DataArray,DataSize):
        """Load varables into class
        example:
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        if(Keyword in SupportKeyWords):#KeyWords Check
            assert len(DataArray)==DataSize,'\n     [Error-%s] Incompatible data size! %d-%d' %(Keyword,len(DataArray),DataSize)
            KeywordID=SupportKeyWords.index(Keyword)
            print('     [%s] '%(Keyword),end='')
            self.SpatialDatas[Keyword]=np.array(DataArray,dtype=KeyWordsDatatypes[KeywordID])
        else:
            if(self.SkipedKeywords==0):print()
            print('     [Warnning] Unsupport keywords[%s]' % (Keyword))
            self.SkipedKeywords+=1


    def read_IncludeFile(self,filename_include,NumData):
        """Read Include data file
        this data file just a series of values
        e.g. 0.2 0.3 12.23 ....
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Aug. 2018
        """

        f=open(filename_include)
        contents=f.read()
        block_dataset=contents.strip().split() #Sepeart input file by slash /
        block_dataset=np.array(block_dataset,dtype=float)
        if(len(block_dataset)!=NumData):
            print('Data size %s is not equal to defined block dimension (NX*NY*NZ) %s'%(len(block_dataset),NumData))
        return block_dataset

    
    def field_cutter(self,nx_range=(0,-1),ny_range=(0,-1),nz_range=(0,-1)):
        """Extract the subset of a domain
        
        Arguments
        ---------
        nx_range    -- The specifc grid range in x for the subset 
        nx_range    -- The specifc grid range in y for the subset 
        nz_range    -- The specifc grid range in z for the subset 
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Feb. 2018

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Feb. 2018
        """

        #If no nx,ny,nz range are defined, all perm will be updated       
        if(nx_range[1]==-1):
            nx_range[1] = self.NX
        if(ny_range[1] == -1):
            ny_range[1] = self.NY
        if(nz_range[1]==-1):
            nz_range[1] = self.NZ

        NX_new=nx_range[1]
        NY_new=ny_range[1]
        NZ_new=nz_range[1]
        N_new=NX_new*NY_new*NZ_new

        TOPS_new=np.zeros(N_new)
        DX_new=np.zeros(N_new)
        DY_new=np.zeros(N_new)
        DZ_new=np.zeros(N_new)

        PERMX_new=np.zeros(N_new)
        PERMY_new=np.zeros(N_new)
        PERMZ_new=np.zeros(N_new)
        PORO_new=np.zeros(N_new)

        ijk_new=0
        for k in range(self.NZ):
            for j in range(self.NY):
                for i in range(self.NX):
                    if(i>=nx_range[0] and i<=nx_range[1]):
                        if(j>=ny_range[0] and j<=ny_range[1]):
                            if(k>=nz_range[0] and k<=nz_range[1]):
                                ijk = getIJK(i, j, k, self.NX, self.NY, self.NZ)
                                DX_new[ijk_new]=self.DX[ijk]
                                DY_new[ijk_new]=self.DY[ijk]
                                DZ_new[ijk_new]=self.DZ[ijk]
                                TOPS_new[ijk_new]=self.TOPS[ijk]

                                PERMX_new[ijk_new]=self.SpatialDatas["PERMX"][ijk]
                                PERMY_new[ijk_new]=self.SpatialDatas["PERMY"][ijk]
                                PERMZ_new[ijk_new]=self.SpatialDatas["PERMZ"][ijk]
                                PORO_new[ijk_new]=self.SpatialDatas["PORO"][ijk]
                                ijk_new=ijk_new+1
        
        NewGrid=GRDECL_Parser('',NX_new,NY_new,NZ_new)
        NewGrid.GRID_type='Cartesian' #Currently only support CartGrid
        NewGrid.DX=DX_new
        NewGrid.DY=DY_new
        NewGrid.DZ=DZ_new
        NewGrid.TOPS=self.TOPS
        NewGrid.SpatialDatas["PERMX"]=PERMX_new
        NewGrid.SpatialDatas["PERMY"]=PERMY_new
        NewGrid.SpatialDatas["PERMZ"]=PERMZ_new
        NewGrid.SpatialDatas["PORO"]=PORO_new
        NewGrid.print_info()
        return NewGrid


    def print_info(self):
        print("     Grid Type=%s Grid" %(self.GRID_type))
        print("     Grid Dimension(NX,NY,NZ): (%s x %s x %s)"%(self.NX,self.NY,self.NZ))
        print("     NumOfGrids=%s"%(self.N))

    ######[DataInterperator]######
    def getPillar(self,Pid):
        """Get a pillar line from COORD
        Pillar is the vertical cell edge line (top point-bottm point)
        
        IndexMap of COORD
        [Row1] xtop ytop ztop xbottom ybottom zbottom
        [Row2] xtop ytop ztop xbottom ybottom zbottom
        ....
        Row follows an order of X->Y->Z

        Arguments
        ---------
        Pid -- Pillar index in [COORD]

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        id_top=[6*Pid+0,6*Pid+1,6*Pid+2]
        id_bottom=[6*Pid+3,6*Pid+4,6*Pid+5]
        TopPoint=np.array([self.COORD[i] for i in id_top])
        BottomPoint=np.array([self.COORD[i] for i in id_bottom])
        return [TopPoint,BottomPoint]

    def getCellPillars(self,i,j):
        """Obtain the four pillars (p0,p1,p2,p3) of a corner point cell
        The index of pillar
        
        3x3x1 system (2D X-Y plane)
        12--- 13  --- 14  ---15
        |      |       |      |  <- Cell 6,7,8
        8 ---  9  --- 10  ---11
        |      |       |      |  <- Cell 3,4,5
        4 ---  5  ---  6  --- 7
        |      |       |      |  <- Cell 0,1,2
        0 ---  1 ---   2 ---  3
        
        
        The pillars index for a grid follows below ordering (XY Plane)
        pil2   pil3
        *------*
        |      |
        |      |
        *------*
        pil0   pil1

        0   12  3
        1. neighboring cell share one common edge index
        2. ZCORN follows the same order for a cell
        3. Try a 3x3x1 grid system using mrst

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        nx,ny=self.NX+1,self.NY+1
        pil0_id,pil1_id=getIJK(i,j,0,nx,ny,0),getIJK(i+1,j,0,nx,ny,0)
        pil2_id,pil3_id=getIJK(i,j+1,0,nx,ny,0),getIJK(i+1,j+1,0,nx,ny,0)

        return [self.getPillar(pil0_id),self.getPillar(pil1_id),self.getPillar(pil2_id),self.getPillar(pil3_id)]     

    def getCellZ(self,i,j,k):
        """Get the Z coords for a cell
        
        Follow getCornerPointCellIdx convention:
        Z, [0,1,2,3,4,5,6,7]

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        CellIds=self.getCornerPointCellIdx(i,j,k)
        return [self.ZCORN[i] for i in CellIds]

    def getCellFaceZ(self,i,j,k,Face='X-,X+,Y-,Y+'):
        """Get the Z coords for a cell
        
         6----7
         -   -   <-Bottom Face
        4----5
          2----3
         -    -  <-Top Face
        0----1   

        Follow getCornerPointCellIdx convention:
        X-, [0,2,4,6]
        X+, [1,3,5,7]
        Y-, [0,1,4,5]
        Y+, [2,3,6,7]

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        CellIds=self.getCornerPointCellIdx(i,j,k)
        if(Face=="X-"): FaceIds=[CellIds[0],CellIds[2],CellIds[4],CellIds[6]]
        if(Face=="X+"): FaceIds=[CellIds[1],CellIds[3],CellIds[5],CellIds[7]]
        if(Face=="Y-"): FaceIds=[CellIds[0],CellIds[1],CellIds[4],CellIds[5]]
        if(Face=="Y+"): FaceIds=[CellIds[2],CellIds[3],CellIds[6],CellIds[7]]

        return [self.ZCORN[i] for i in FaceIds]

    def getCellCoords(self,i,j,k):
        """Get XYZ coords for eight node of a cell

        6----7
         -   -   <-Bottom Face
        4----5
          2----3
         -    -  <-Top Face
        0----1   

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        XYZ=[]
        Pillars=self.getCellPillars(i,j)
        Zs=self.getCellZ(i,j,k)
        for pi in range(8): # Loop 8 point for each cell
            Pillar_i=pi%4
            XYZ.append(self.interpPtsOnPillar(Zs[pi],Pillars[Pillar_i]))
        return XYZ
    
    def getCell1Coord(self,i,j,k,nodeID):
        """Get XYZ coords for one node of 8 nodes for a cell

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        Pillars=self.getCellPillars(i,j)
        Zs=self.getCellZ(i,j,k)
        Pillar_i=nodeID%4
        return self.interpPtsOnPillar(Zs[nodeID],Pillars[Pillar_i])


    def getCornerPointCellIdx(self,i,j,k):
        """Obtain the eight coords index for a cell

        3x3x1 system (2D X-Y plane)
        30---31,32---33,34---35
        |      |       |      |  <- Cell 6,7,8
        24---25,26---27,28---29
        18---19,20---21,22---23
        |      |       |      |  <- Cell 3,4,5
        12---13,14---15,16---17
        6 --- 7,8 --- 9,10---11
        |      |       |      |  <- Cell 0,1,2
        0 --- 1,2 --- 3,4 --- 5

        Node order convention for a 3D cell
         6----7
         -   -   <-Bottom Face
        4----5
          2----3
         -    -  <-Top Face
        0----1

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        nx,ny,nz=2*self.NX,2*self.NY,2*self.NZ

        p1_id,p2_id=getIJK(2*i,2*j,2*k,nx,ny,nz),getIJK(2*i+1,2*j,2*k,nx,ny,nz)
        p3_id,p4_id=getIJK(2*i,2*j+1,2*k,nx,ny,nz),getIJK(2*i+1,2*j+1,2*k,nx,ny,nz)

        p5_id,p6_id=getIJK(2*i,2*j,2*k+1,nx,ny,nz),getIJK(2*i+1,2*j,2*k+1,nx,ny,nz)
        p7_id,p8_id=getIJK(2*i,2*j+1,2*k+1,nx,ny,nz),getIJK(2*i+1,2*j+1,2*k+1,nx,ny,nz)

        #print(p1_id,p2_id,p3_id,p4_id)#Top Layer
        #print(p5_id,p6_id,p7_id,p8_id)#Bottom Layer

        return p1_id,p2_id,p3_id,p4_id,p5_id,p6_id,p7_id,p8_id

    def interpPtsOnPillar(self,z,Pillar):
        """Obtain the eight coords for a cell
           X,Y coords has to be interpolated from Z
        xy1=xy0+k*z

        Pillar=[(x0 y0 z0),(x1 y1 z1)]
        (x,y,z) is somewhere between (x0 y0 z0) and (x1 y1 z1)

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        if(abs(Pillar[1][2]-Pillar[0][2])>1e-8):
            k=(z-Pillar[0][2])/(Pillar[1][2]-Pillar[0][2])
        else:#Degenrated cell
            k=0.0

        x=Pillar[0][0]+k*(Pillar[1][0]-Pillar[0][0])
        y=Pillar[0][1]+k*(Pillar[1][1]-Pillar[0][1])
        
        return np.array([x,y,z])
    
    def detectFaceFault(self,Z_ijk,Z_neigh):
        """#* Check Fault type for a face (X-,X+,Y-,Y+)
        TypeID
        -1  --Fully connected
        0  --Sealing
        >0--Partially Fault, fault gap value, e.g. 1.5..

        Z_ijk face node order
        p2     p3
        *------*
        |      |
        |      |
        *------*
        p0     p1

        Arguments
        ---------
        Z_ijk   -- z value of a cell face,e.g. [0.5,0.2,0.6,0.7]
        Z_neigh -- z value of a neighbor cell face

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """

        '''
        overlap_p02=overlap(Z_ijk[0],Z_ijk[2],Z_neigh[0],Z_neigh[2])
        overlap_p13=overlap(Z_ijk[1],Z_ijk[3],Z_neigh[1],Z_neigh[3])

        length_p02=Z_ijk[2]-Z_ijk[0]
        length_p13=Z_ijk[3]-Z_ijk[1]

        gap_p02=length_p02-overlap_p02
        gap_p13=length_p13-overlap_p13

        print('Overlap',overlap_p02,overlap_p13)
        print('Gap',gap_p02,gap_p13)
       

        if(abs(gap_p02)+abs(gap_p13)<1e-10): #Fully connected
           
            return -1.0
        elif(abs(overlap_p02)+abs(overlap_p13)<1e-10): #Sealing fault
            return 0.0
        else:#Partially connected
            return 0.5
        '''
        #Simple method
        diffVec=np.array(Z_ijk)-np.array(Z_neigh)
        if(sum(abs(diffVec))>0): 
            return 0.5
        else:
            return -1


    def isBoundaryCell(self,Cell=[0,0,0],Dim='3D'):
        ''' Check the a given cell is boundary cell or not

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        '''
        count=0
        face=[]
        #Boundary Point
        if(Cell[0]==0):
            count+=1
            face.append('X-')
        if(Cell[0]==self.NX-1):
            count+=1
            face.append('X+')
        if(Cell[1]==0):
            count+=1
            face.append('Y-')
        if(Cell[1]==self.NY-1):
            count+=1
            face.append('Y+')

        if(Dim=="3D"):
            if(Cell[2]==0):
                count+=1
                face.append('Z-')
            if(Cell[2]==self.NZ-1):
                count+=1
                face.append('Z-')

        return count,face

    def findCellFault(self,Cell=[0,0,0]):
        ''' Check the fault for 4 faces of a cell [X-,X+,Y-,Y+] 2D
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        '''
        i,j,k=Cell
        Faces=['X-','X+','Y-','Y+']
        Fault=[False,False,False,False]
        
        FaultMarker=-1
        Z_ijk=self.getCellFaceZ(i,j,k,'X-')
        if(i!=0): 
            Z_neigh=self.getCellFaceZ(i-1,j,k,'X+')
            FaultMarker=self.detectFaceFault(Z_ijk,Z_neigh)
            #print(Z_ijk,Z_neigh,FaultMarker)
        if(FaultMarker!=-1.0 or i==0):#This is a fault here
            Fault[0]=True
        
        FaultMarker=-1
        Z_ijk=self.getCellFaceZ(i,j,k,'X+')
        if(i!=self.NX-1): 
            Z_neigh=self.getCellFaceZ(i+1,j,k,'X-')
            FaultMarker=self.detectFaceFault(Z_ijk,Z_neigh)
        if(FaultMarker!=-1.0 or i==self.NX-1):#This is a fault here
            Fault[1]=True

        FaultMarker=-1
        Z_ijk=self.getCellFaceZ(i,j,k,'Y-')
        if(j!=0):
            Z_neigh=self.getCellFaceZ(i,j-1,k,'Y+')
            FaultMarker=self.detectFaceFault(Z_ijk,Z_neigh)
        if(FaultMarker!=-1.0 or j==0):#This is a fault here
            Fault[2]=True
        
        FaultMarker=-1
        Z_ijk=self.getCellFaceZ(i,j,k,'Y+')
        if(j!=self.NY-1): 
            Z_neigh=self.getCellFaceZ(i,j+1,k,'Y-')
            FaultMarker=self.detectFaceFault(Z_ijk,Z_neigh)
            #print(Z_ijk,Z_neigh,FaultMarker)
        if(FaultMarker!=-1.0 or j==self.NY-1):#This is a fault here
            Fault[3]=True

        return Fault

#############################################
#
#  Auxiliary function
#
#############################################

def parseDataArray(DataArray):
        """Parse special dataArray format in GRDECL 
        example:
            5*3.0=[3.0 3.0 3.0 3.0 3.0]
            1.0 2*3.0 5.0=[1.0 3.0 3.0 5.0]
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """

        data=[]
        error_count=0
        for value in DataArray:
            if(is_number(value)==2):
                num,val=value.split('*')
                for i in range(int(num)): data.append(val)
            elif(is_number(value)==1):
                data.append(value)
            else:
                error_count+=1
        
        if(error_count>0):
            print(DataArray)
        
        assert error_count==0, '[Error] Can not find any numeric value!'
        
        

        return data

def KeyWordReader(fname,varname, datatype=float):
    #Simple reader to read a file with a input keyword name

    f=open(fname)
    print('[Input] Reading ECLIPSE/PETREL file \"%s\" ....'%(fname))
    contents=f.read()
    f.close()

    contents=RemoveCommentLines(contents,commenter='--')
    contents_in_block=contents.strip().split('/') #Sepeart input file by slash /
    contents_in_block = [x for x in contents_in_block if x]#Remove empty block at the end
    NumKeywords=len(contents_in_block)

    for i,block in enumerate(contents_in_block):#Keyword, Block-wise
        #Clean the data where no spliter \ provided
        block=scanKeyword(block)

        blockData_raw=block.strip().split()
        Keyword,DataArray=blockData_raw[0],blockData_raw[1:]

        if(Keyword==varname):
            print("     Reading Keywords [%s] " %(Keyword))
            DataArray=parseDataArray(DataArray)
            return np.array(DataArray,dtype=datatype)

    print('     [Warnning] Can not find keywords on file[%s]',varname)

    return None


def RemoveCommentLines(data,commenter='--'):
    #Remove comment and empty lines
    data_lines=data.strip().split('\n')
    newdata=[]
    for line in data_lines:
        if line.startswith(commenter) or not line.strip():
            # skip comments and blank lines
            continue   
        newdata.append(line)
    return '\n'.join(newdata)

def scanKeyword(data):
    #scan and find the keyword
    #e.g. ['INIT','DX','2500*30.0'] -> ['DX','2500*30.0']
    for key in SupportKeyWords:
        if (key in data) and (data.find(key)!=0):
            return data[data.find(key):-1]
    return data

def is_number(s):
    #Determine a string is a number or not
    #Used in [read_GRDECL] [getBlkdata]
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    
    try: #Special format N*val= [val val val ....]
        num, val = s.split('*')
        return 2
    except ValueError:
        pass
 
    return False

def getI_J_K(ijk,NX,NY,NZ):
    #Find index [i,j,k] from a flat 3D matrix index [i,j,k]
    i,j,k=0,0,0
    
    i=ijk%NX
    j=((int)(ijk / NX)) % NY
    k=(int)(ijk / (NX*NY))
    
    return i,j,k

def getIJK(i,j,k,NX,NY,NZ):
    #Convert index [i,j,k] to a flat 3D matrix index [ijk]
    return i + (NX)*(j + k*(NY))

def overlap(min1, max1, min2, max2):
    #Math: overlap of 1D segments
    #Reference: https://stackoverflow.com/questions/16691524/calculating-the-overlap-distance-of-two-1d-line-segments?rq=1
    return max(0.0, min(max1, max2) - max(min1, min2))
