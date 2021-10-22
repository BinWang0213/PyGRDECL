import numpy as np


def getIJK(i,j,k,NX,NY,NZ):
    #Convert index [i,j,k] to a flat 3D matrix index [ijk]
    return i + (NX)*(j + k*(NY))

def Cartesian2UnstructGrid(DX,DY,DZ,TOPS,NX,NY,NZ):
    #Convert self.DX to coordiantes-X of unstructure gridblocks
    #Used in [write_VTU]
    debug=0
    
    coordX,coordY,coordZ=np.zeros((2*NX,2*NY,2*NZ)),np.zeros((2*NX,2*NY,2*NZ)),np.zeros((2*NX,2*NY,2*NZ))
    for k in range(2*NZ):
        for j in range(2*NY):
            for i in range(2*NX):
                #print(DX[i][j])
                i_GB,j_GB,k_GB=(int)(i/2),(int)(j/2),(int)(k/2)

                #CoordX
                if(i==0):
                    coordX[0][j][k]=0
                    if(debug):print('First Node',(0),'-')
                if(i>0 and i%2==0):
                    ijk_GB=getIJK(i_GB-1,j_GB,k_GB,NX,NY,NZ)
                    coordX[i-1][j][k]=DX[ijk_GB] #odd
                    coordX[i][j][k]=DX[ijk_GB] #even
                    if(debug):print('pair',(i-1,i),(i_GB-1,j_GB))
                    if(i_GB>1):
                        coordX[i-1][j][k]+=coordX[i-1-1][j][k]
                        coordX[i][j][k]=coordX[i-1][j][k]
                if(i==2*NX-1):
                    ijk_GB=getIJK(NX-1,j_GB,k_GB,NX,NY,NZ)
                    coordX[2*NX-1][j][k]=DX[ijk_GB]+coordX[2*NX-1-1][j][k]
                    if(debug):print('Last Node',(2*NX-1),(NX-1,j_GB))
                
                #CoordY
                if(j==0):
                    coordY[i][0][k]=0
                    if(debug and i==0 and k==0):print('First Node',(0),'-')
                if(j>0 and j%2==0):
                    ijk_GB=getIJK(i_GB,j_GB-1,k_GB,NX,NY,NZ)
                    coordY[i][j-1][k]=DY[ijk_GB] #odd
                    coordY[i][j][k]=DY[ijk_GB] #even
                    if(debug and i==0 and k==0):print('pair',(j-1,j),(i_GB,j_GB-1))
                    if(j_GB>1):
                        coordY[i][j-1][k]+=coordY[i][j-1-1][k]
                        coordY[i][j][k]=coordY[i][j-1][k]
                if(j==2*NY-1):
                    ijk_GB=getIJK(i_GB,NY-1,k_GB,NX,NY,NZ)
                    coordY[i][2*NY-1][k]=DY[ijk_GB]+coordY[i][2*NY-1-1][k]
                    if(debug and i==0 and k==0):print('Last Node',(2*NY-1),(i_GB,NY-1))
                
                #CoordZ
                ijk_GB=getIJK(i_GB,j_GB,k_GB,NX,NY,NZ)
                
                if(k==0):
                    coordZ[i][j][0]=TOPS[ijk_GB]
                    if(debug and i==0 and j==0):print('First Node',(0),'-')
                if(k>0 and k%2==0):
                    ijk_GB=getIJK(i_GB,j_GB,k_GB-1,NX,NY,NZ)
                    coordZ[i][j][k-1]=DZ[ijk_GB] #odd
                    coordZ[i][j][k]=DZ[ijk_GB] #even
                    if(debug and i==0 and j==0):print('pair',(k-1,k),(i_GB,k_GB-1))
                    if(k_GB>1):
                        coordZ[i][j][k-1]+=coordZ[i][j][k-1-1]
                        coordZ[i][j][k]=coordZ[i][j][k-1]
                if(k==2*NZ-1):
                    ijk_GB=getIJK(i_GB,j_GB,NZ-1,NX,NY,NZ)
                    coordZ[i][j][2*NZ-1]=DZ[ijk_GB]+coordZ[i][j][2*NZ-1-1]
                    if(debug and i==0 and j==0):print('Last Node',(2*NZ-1),(i_GB,NZ-1))

    return coordX,coordY,coordZ

#Create random permeability field
def logNormLayers(gridDims,AvgLayerPerm,poro_const=0.05):
    import scipy.ndimage as ndimage
    import matplotlib.pyplot as plt
    
    assert len(AvgLayerPerm)==gridDims[2], "[Error] AvgLayerPerm size is not equal to model layers!\n"
    
    NumLayers=gridDims[2]
    K=[]
    for li in range(NumLayers):
        #Genetrate random field using gaussian smmoth, see logNormLayers.M in MRST
        raw_data = np.random.randn(gridDims[0],gridDims[1])
        raw_data -= 0.05*np.random.randn(1)[0]
        k = ndimage.gaussian_filter(raw_data, 3.5)
        k=np.exp( 2 + 2*k)
        
        K.append(AvgLayerPerm[li]*k/np.mean(k))
        print("Layer%d Perm Avg=%lf min=%lf max=%lf" %(li+1,np.average(K[li]),np.min(K[li]),np.max(K[li])))
    
    #Plot the permeability
    NumRows=int(NumLayers/3)+1
    fig, axs = plt.subplots(NumRows,3, figsize=(4*3, NumRows*3), facecolor='w', edgecolor='k')
    #fig.subplots_adjust(hspace = .01, wspace=.001)

    axs = axs.ravel()
    for li in range(NumLayers):
        logperm=np.log10(K[li])
        #logperm=K[li]
        im=axs[li].imshow(logperm, interpolation='nearest',origin='lower',cmap='jet',
                          #vmin=np.min(K),vmax=np.max(K))
                          vmin=np.min(np.log10(K)),vmax=np.max(np.log10(K)))
        axs[li].set(title='Random Log(perm) Layer%d'%(li+1), xticks=[], yticks=[])

    for ax in axs:
        ax.set(xticks=[], yticks=[])

    
    #fig.tight_layout()
    cbar = fig.colorbar(im, ax=axs.ravel().tolist(), shrink=0.95)
    plt.show()
    
    #Genetrate random porosity field using phi=0.25*K^0.11
    phi=[]
    for li in range(NumLayers):
        phi.append(poro_const*np.log(K[li]))
        print("Layer%d Poro Avg=%lf min=%lf max=%lf" %(li+1,np.average(phi[li]),np.min(phi[li]),np.max(phi[li])))

    #Plot the Porosity
    #NumRows=max([int(NumLayers/3),1])
    fig, axs = plt.subplots(NumRows,3, figsize=(4*3, NumRows*3), facecolor='w', edgecolor='k')
    #fig.subplots_adjust(hspace = .01, wspace=.001)

    axs = axs.ravel()
    for li in range(NumLayers):
        im=axs[li].imshow(phi[li], interpolation='nearest',origin='lower',cmap='jet',
                          vmin=np.min(phi),vmax=np.max(phi))
        axs[li].set(title='Random Porosity Layer%d'%(li+1), xticks=[], yticks=[])
    
    for ax in axs:
        ax.set(xticks=[], yticks=[])
    
    
    #fig.tight_layout()
    cbar = fig.colorbar(im, ax=axs.ravel().tolist(), shrink=0.95)
    plt.show()
    
    #Flatten K and phi into x,y,z GRDECL convention order
    return np.array([i for k in K for i in k.flatten()]),np.array([i for fi in phi for i in fi.flatten()])