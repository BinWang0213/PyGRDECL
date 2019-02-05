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

