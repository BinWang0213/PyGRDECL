#########################################################################
#       (C) 2021 Department of Numerical Geology,                       #
#       Universite de lorraine, CNRS, OTELO, Nancy, France.             #
#                                                                       #
# This code is released under the terms of the BSD license, and thus    #
# free for commercial and research use. Feel free to use the code into  #
# your own project with a PROPER REFERENCE.                             #
#                                                                       #
# Permeability Tools for PyGRDECL Code (Bin Wang)                       #
# Author: Zakari Mustapha                                               #
# Email: mustapha.zakari@univ-lorraine.fr                               #
#########################################################################

import numpy as np
import math

# MZ:Generate Fine scale random K, scalar values
def gen_rand_K(Grid):
#   1 mD = 9.871e-16;
    mean=0;
    std=1;
    size=(Grid.NX,Grid.NY,Grid.NZ)
    K=np.empty(size,dtype='float')
    K=np.random.normal(mean,std,size=size)
    K=np.exp(K*math.log(10)) # 10^K
    return K.reshape(np.product(size), order='F')

# MZ: Plot histogram at logscale
import matplotlib.pyplot as plt
def plot_hist(Kvect,varname=""):
    hist, bins = np.histogram(Kvect, bins=30);
    # Lognormal distribution
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    
    fig, axs = plt.subplots(1,1, figsize=(10, 5), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace = .5, wspace=.001)

    axs.set_title(varname)     
    axs.hist(Kvect, bins=logbins, edgecolor='black', linewidth=1.2)
    plt.xscale('log')
    plt.show()
    
# MZ::Add TPFA Pressure calculation
def compute_TPFA_Pressure(Grid, K):
    import scipy.sparse as SP
    from scipy.sparse.linalg import spsolve
    Nx=Grid.NX ; Ny=Grid.NY ; Nz=Grid.NZ
    N = Nx*Ny*Nz
    Dx = Grid.DX;    Dy = Grid.DY;    Dz = Grid.DZ

    q = 0*np.ones_like(Dx)
    Dx = Dx.reshape((Nx, Ny, Nz), order='F')
    Dy = Dy.reshape((Nx, Ny, Nz), order='F')
    Dz = Dz.reshape((Nx, Ny, Nz), order='F')

    # Inverse K values stored in L
    K=np.ones_like(K)
    invK=K**(-1);#np.ones_like(K);
    DxinvK=Dx*invK
    DyinvK=Dy*invK
    DzinvK=Dz*invK

    Ord='F'

    # Initialize transmissibilities
    # Interface area and distance to it
    Areax=Dy*Dz;
    Areay=Dx*Dz;
    Areaz=Dy*Dx;
    TX=np.zeros((Nx+1,Ny,Nz));
    TY=np.zeros((Nx,Ny+1,Nz));
    TZ=np.zeros((Nx,Ny,Nz+1));

    TX[1:Nx,:,:]=2*(Areax[0:Nx-1,:,:]+ Areax[1:Nx,:,:])/\
                  (DxinvK[0:Nx-1,:,:]+DxinvK[1:Nx,:,:]) #no contrib in xmin xmax
    TY[:,1:Ny,:]=2*(Areay[:,0:Ny-1,:]+ Areay[:,1:Ny,:])/\
                  (DyinvK[:,0:Ny-1,:]+DyinvK[:,1:Ny,:]) #no contrib in ymin ymax
    TZ[:,:,1:Nz]=2*(Areaz[:,:,0:Nz-1]+ Areaz[:,:,1:Nz])/\
                  (DzinvK[:,:,0:Nz-1]+DzinvK[:,:,1:Nz]) #no contrib in zmin zmax

    # Assembling
    x1=TX[0:Nx,:,:].reshape((N), order=Ord);  x2=TX[1:Nx+1,:,:].reshape((N), order=Ord)
    y1=TY[:,0:Ny,:].reshape((N), order=Ord);  y2=TY[:,1:Ny+1,:].reshape((N), order=Ord)
    z1=TZ[:,:,0:Nz].reshape((N), order=Ord);  z2=TZ[:,:,1:Nz+1].reshape((N), order=Ord)

    diagonals = np.zeros((7, N))
    diagonals[0,0:N] = -z2;
    diagonals[1,0:N] = -y2;
    diagonals[2,0:N] = -x2;
    diagonals[3,0:N] = x1+x2+y1+y2+z1+z2;
    diagonals[4,0:N] = -x1;
    diagonals[5,0:N] = -y1;
    diagonals[6,0:N] = -z1;

    dec=[-Ny*Nx,-Nx,-1,0,1,Nx,Ny*Nx]
    A = SP.spdiags(diagonals, dec, N, N, format='lil')
    ind=0
    # for ind in bdry.ind_L_:
    A[ind,:]=0
    A[ind,ind]=1;
    q[ind]=-1

    ind=N-1
    # for ind in bdry.ind_L_:
    A[ind,:]=0
    A[ind,ind]=1;
    q[ind]=1
    # for ind in bdry.ind_R_:
    #     A[ind,:]=0
    #     A[ind,ind]=1
    #     q[ind]=bdry.val_R_


    A=A.tocsr()
    p = spsolve(A, -q)
    return p