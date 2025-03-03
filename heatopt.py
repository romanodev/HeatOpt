import numpy as np
from itertools import product
import nlopt
from functools import partial
from scipy.signal import convolve2d as convolve_scipy
import scipy
import scipy.sparse as sp
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import warnings
import matplotlib.pylab as plt
from stl import mesh,Mode
from  io import BytesIO
import matplotlib as mpl
import scipy.sparse.linalg as spla
import time
import numpy.ma as ma
import json

#mpl.rcParams['toolbar'] = 'None'

from scipy.sparse import (spdiags, SparseEfficiencyWarning, csc_matrix,
    csr_matrix, isspmatrix, dok_matrix, lil_matrix, bsr_matrix)
warnings.simplefilter('ignore',SparseEfficiencyWarning)




def create_stl(x):

 x = np.array(x)   

 N = len(x)
 Ns = int(np.sqrt(N))
 x =x.reshape((Ns,Ns))

 #Make it periodic--
 x = np.pad(x,Ns,mode='wrap')

 N2 = 3*Ns

 vertices = []
 faces    = []

 for i, j in product(range(N2+1),range(N2+1)):
     vertices.append([i,j,0])
 for i, j in product(range(N2+1),range(N2+1)):
     vertices.append([i,j,1])

 def add_pixel(i,j):
    v  = [] 
    v.append(i*(N2+1)+j)
    v.append((i+1)*(N2+1)+j)
    v.append((i+1)*(N2+1)+j+1)
    v.append(i*(N2+1)+j+1)
    v.append(i*(N2+1)+j + (N2+1)*(N2+1))
    v.append((i+1)*(N2+1)+j + (N2+1)*(N2+1))
    v.append((i+1)*(N2+1)+j+1 + (N2+1)*(N2+1))
    v.append(i*(N2+1)+j+1 + (N2+1)*(N2+1))
   

    faces.append([v[0],v[3],v[1]])
    faces.append([v[1],v[3],v[2]])
    faces.append([v[0],v[4],v[7]])
    faces.append([v[0],v[7],v[3]])
    faces.append([v[4],v[5],v[6]])
    faces.append([v[4],v[6],v[7]])
    faces.append([v[5],v[1],v[2]])
    faces.append([v[5],v[2],v[6]])
    faces.append([v[2],v[3],v[6]])
    faces.append([v[3],v[7],v[6]])
    faces.append([v[0],v[1],v[5]])
    faces.append([v[0],v[5],v[4]])


 for i,j in zip(*x.nonzero()):
  add_pixel(i,j)    

 faces    = np.array(faces)
 vertices = np.array(vertices)

 T = 4
 vertices[:,2] *=T


 # Create the mesh
 cube = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
 for i, f in enumerate(faces):
    for j in range(3):
        cube.vectors[i][j] = vertices[f[j],:]

 output = BytesIO()
 cube.save('TopOpt',fh=output,mode=Mode.ASCII)

 # Write the mesh to file "cube.stl"
 return output.getvalue().decode("utf-8") 


def get_fig(x,grid,J):

  x = np.where(x>0.5,1,0).reshape((grid,grid))
  
  data = np.linalg.norm(J,axis=(0,1)).reshape((grid,grid))
 
  mask = x==0
  mask = np.pad(mask.T,grid,mode='wrap')

  data = np.pad(data.T,grid,mode='wrap')

  data =  ma.array(data,mask=mask)
  
  #cmap = ListedColormap(["#cadbde","#404040"])
  cmap = 'viridis'

  if type(x) == list: x = np.array(x)  
  fig  = plt.figure(facecolor='#404040')
  ax   = fig.add_subplot(111)
  #ax.imshow(1-np.pad(x.reshape((grid,grid)).T,grid,mode='wrap'),vmin=0,vmax=1,cmap=cmap,animated=True)

  data -= data.min()
  data /=data.max()
  im = ax.imshow(data,cmap=cmap,animated=True,vmin=0,vmax=1)
  ax.plot([grid,grid,2*grid,2*grid,grid],[grid,2*grid,2*grid,grid,grid],color='white',ls='--',lw=2)
  ax.axis('off')

  # Add a colorbar at the bottom
  cbar = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)
  cbar.set_label('Magnitude of Heat Flux (Normalized)',color='white')
  cbar.ax.xaxis.set_tick_params(color='white')
  cbar.ax.tick_params(labelcolor='white')
  cbar.outline.set_edgecolor('none')
  

  fig.tight_layout(pad=0,h_pad=0,w_pad=0)

  return fig


def get_grid(N):


    def maps(i,j):

        return (i%N)*N + j%N
    #Compute indices for adiacent elements
    i_mat = []
    j_mat = []
    normals = []
    k  = 0
    for i in range(N):
     for j in range(N):

         k1,k2 = maps(i,j),maps(i+1,j)

         i_mat.append(k1)    
         j_mat.append(k2)    
         normals.append([1,0])
         k +=1

         k2 = maps(i,j+1)
         i_mat.append(k1)    
         j_mat.append(k2)    
         normals.append([0,-1])
         k +=1

         k2 = maps(i,j-1)
         i_mat.append(k1)    
         j_mat.append(k2)    
         normals.append([0,1])
         k +=1

         k2 = maps(i-1,j)
         i_mat.append(k1)    
         j_mat.append(k2)    
         normals.append([-1,0])
         k +=1
         
    #These relationships are found heuristically      

    kr = 4*N*N-4*N + 4*np.arange(N)  
    ku = 2+4*N*np.arange(N)
    kd = 4*N-3 +4*N*np.arange(N)
    kl = 3+4*np.arange(N)

    #calculate centroids--
    x = np.linspace(-1/2+1/(2*N),1/2-1/(2*N),N)
    centroids = np.array([ [i,j] for i in x for j in x[::-1]])
    #--------------------
    ind_extremes = [[kl,kr],[ku,kd]]

    return {'i':np.array(i_mat),'j':np.array(j_mat),'ind_extremes':np.array(ind_extremes,dtype=int),'normals':np.array(normals),'centroids':centroids}





def conic_filter_2D(centroids,R):

    N = np.array(np.sqrt(len(centroids)),int)

    tmp = np.sqrt(np.power(centroids[:,0],2) + np.power(centroids[:,1],2))/R

    conic_kernel = np.where(tmp<1,1-tmp,0).reshape((N,N))

    if not (np.sum(conic_kernel) == 0):
      conic_kernel /= np.sum(conic_kernel)

    return conic_kernel

def get_mapping(grid,R):

   N = int(grid**2.0)
   shift = int(grid/2)-1
   eta = 0.5
   aux     = get_grid(grid)
   kernel  = conic_filter_2D(aux['centroids'],R)

   gt = np.zeros((grid,grid,grid,grid))
   for i1 in np.arange(-shift,grid-shift):
     for j1 in np.arange(-shift,grid-shift):
      gt[i1,j1] = kernel[(i1+shift)%grid-np.arange(grid)[:,np.newaxis],(j1+shift)%grid-np.arange(grid)[np.newaxis,:]]

   gt2 = gt.reshape((N,N))

   def mapping(x,beta):

    xt = convolve_scipy(x.reshape((grid,grid)),kernel,mode='same',boundary='wrap').flatten()

    xp = (np.tanh(beta*eta) + np.tanh(beta*(xt-eta)))/(np.tanh(beta*eta) + np.tanh(beta*(1-eta)))

    proj_numpy_jac = beta*(1-np.power(np.tanh(beta*(xt-eta)),2))/(np.tanh(beta*eta) + np.tanh(beta*(1-eta)))
   
    gp = np.einsum('ij,i->ij',gt2,proj_numpy_jac)

    return xp,gp

   return mapping


def fourier(**options):
    """Fourier Solver"""

    #Directions---
    #phi = [0,np.pi/2]
    #directions = np.array([np.cos(phi),np.sin(phi)]).T
    #Ainv = np.array([[1,0],[0,1]])

    #This is needed in order to get the whole tensor
    A = np.array([[1,0,0],[0,1,0],[0.5,0.5,1]])
    Ainv = np.linalg.inv(A)
    phi = [0,np.pi/2,np.pi/4]
    directions = np.array([np.cos(phi),np.sin(phi)]).T



    direct = False
    #directions =options['directions']
    n_dir,dim = np.array(directions).shape
    grid        = options['grid']
    if dim == 3:
     N           = int(grid**3)
     factor      = 1/grid
    else: 
     N           = int(grid**2)
     factor      = 1

    aux         = get_grid(grid)
    i_mat       = aux['i']
    j_mat       = aux['j']
    ind_extremes = aux['ind_extremes'] 
    ii = [ind_extremes[i,1]  for i in range(dim)]
    normals = aux['normals']

    #For direct
    row = np.hstack((i_mat,np.arange(N))) 
    col = np.hstack((j_mat,np.arange(N))) 

    k0 = 1e-12;k1 = 1.0

    def func(rho):
    
    
     kappa_map       = k0 + rho*(k1-k0)
     kappa_ij     = 2*kappa_map[i_mat] * kappa_map[j_mat]/(kappa_map[i_mat] + kappa_map[j_mat])

     kappad_ij    = (k1-k0)*0.5*np.power(kappa_ij/kappa_map[i_mat],2)
     kappad_ji    = (k1-k0)*0.5*np.power(kappa_ij/kappa_map[j_mat],2)


     #compute vectorial perturbation
     P_vec   = np.zeros((dim,N))
     for i in range(dim):
      r  = ii[i] 
      np.add.at(P_vec,(i,i_mat[r]), kappa_ij[r])
      np.add.at(P_vec,(i,j_mat[r]),-kappa_ij[r])   



     #Assemble matrix
     tmp = np.zeros_like(rho)
     np.add.at(tmp,i_mat,kappa_ij)
     d   = np.hstack((-kappa_ij,tmp))
     A   = scipy.sparse.csc_matrix((d, (row, col)), shape=(N, N))
     A[12,:] = 0;A[12,12] = 1
     lu = sp.linalg.splu(A)
     #---------------------

     #---------------------------------------
    
    

     P =  np.einsum('di,ic->cd',directions,P_vec)   
     P[12,:] = 0 
     T_mat  = lu.solve(P)
   

     jacobian = np.zeros((n_dir,N))
     np.add.at(jacobian.T,i_mat,np.einsum('s,sd->sd',kappad_ij,np.power(T_mat[i_mat]-T_mat[j_mat],2)))
     tmp = (directions.T**2)[:,None,:] + 2*directions.T[:,None,:]*(T_mat[j_mat[ii]]-T_mat[i_mat[ii]])
     np.add.at(jacobian.T,i_mat[ii],kappad_ij[ii][:,:,None]*tmp)
     np.add.at(jacobian.T,j_mat[ii],kappad_ji[ii][:,:,None]*tmp)

    
     kappa = np.einsum('is,di->d',kappa_ij[ii],directions**2) -  np.einsum('cd,ic,di->d',T_mat,P_vec,directions) 

     
     kappa    = np.matmul(Ainv,kappa)
     #print(kappa)
     jacobian = np.einsum('ji,ik->jk',Ainv,jacobian)


     return kappa,jacobian,T_mat

    
    def compute_flux(rho,T_mat):
      

     kappa_map       = k0 + rho*(k1-k0)
     kappa_ij     = 2*kappa_map[i_mat] * kappa_map[j_mat]/(kappa_map[i_mat] + kappa_map[j_mat])

     J = np.zeros((n_dir,2,N))

     for n,direction in enumerate(directions):
      
      #Flux-----------------------------
      #Bulk
      T_boundary       = (T_mat[i_mat,n]*kappa_map[i_mat] + T_mat[j_mat,n]*kappa_map[j_mat])/(kappa_map[i_mat]+kappa_map[j_mat]) 
      contribution_i   = np.einsum('s,sd->sd',T_boundary*kappa_map[i_mat],normals)
      contribution_j   = np.einsum('s,sd->sd',T_boundary*kappa_map[j_mat],normals)
      np.add.at(J[n,0],i_mat,-contribution_i[:,0])
      np.add.at(J[n,1],i_mat,-contribution_i[:,1])
      np.add.at(J[n,0],j_mat,contribution_j[:,0])
      np.add.at(J[n,1],j_mat,contribution_j[:,1])

     
      #Periodic

      for i in range(dim):
       r = ii[i]

       contribution = direction[i]*np.einsum('s,sd->sd',kappa_ij[r],normals[r])
      
       np.add.at(J[n,0],i_mat[r],-contribution[:,0])
       np.add.at(J[n,1],i_mat[r],-contribution[:,1])
       np.add.at(J[n,0],j_mat[r],-contribution[:,0])
       np.add.at(J[n,1],j_mat[r],-contribution[:,1])



     return J 
    
  #---------------------------------------------------------------



   
    return func,compute_flux


def get_optimizer():

 nb = 14
 mm = 30
 R  = 0.1
 grid    = 40  #Grid resolution
 N = grid**2
 betas = [2**(n+1) for n in range(nb)]
 maxiter = [mm]*nb
 maxiter.append(1)
 betas.append(1e24)

 f,get_flux = fourier(grid=grid)

 mapping = get_mapping(grid,R=R)

 counts = [0]
 def optimize(*kd):
 
  kd = np.array(kd)  
  n_dir = len(kd)
  kappa = np.zeros_like(kd)

  T_mat = np.zeros((N,n_dir))
  def func(x,grad,beta):
   counts[0] +=1   
   x,gp =  mapping(x,beta)
   kappa[:],jacobian,T_mat[:,:] = f(x)
   g = np.linalg.norm(kappa-kd)
   jacobian = np.einsum('ik,kl->il',jacobian,gp)
   
   grad[:] = np.zeros(N)
   for i in range(n_dir):
    grad[:]  += 1/g*((kappa[i]-kd[i])*jacobian[i])

   return g

  opt = nlopt.opt(nlopt.LD_CCSAQ,N)

  x = np.random.rand(N)
  
  #x.dump('x')
  #x = np.load('x',allow_pickle=True)
  #---------------------
  total_count = 0
  for miter,beta in zip(*(maxiter,betas)):

    counts[0] = 0
    opt = nlopt.opt(nlopt.LD_CCSAQ,N)
    opt.set_lower_bounds(np.zeros(N))
    opt.set_upper_bounds(np.ones(N))
    opt.set_min_objective(partial(func,beta=beta))
    opt.set_maxeval(miter)
    #opt.set_xtol_rel(1e-3)
    opt.set_stopval(0.005)
    #opt.set_ftol_abs(0.001)
    x = opt.optimize(x)
    total_count += counts[0]

  x = mapping(x,beta)[0]

  #Prune structure--- 
  a = f(x)
  gradient_norm = np.linalg.norm(a[1],axis=0)
  x[gradient_norm < 1e-7] = 0
  #------------------


  J = get_flux(x,T_mat)
  

  fig = get_fig(x,grid,J)

  #print('total_count')
  #print(total_count)
  return kappa,fig,x,J

 return optimize,grid



if __name__ =='__main__':

 optimizer,grid= get_optimizer()
 a = time.time()
 kappa,fig,x,J = optimizer(0.2,0.2,-0.05)
 print(time.time()-a)

 #Save the json file
 print(kappa)
 data = {'kappa':kappa.tolist(),'x':x.tolist(),'J':J.tolist()}

 with open('example.json','w') as f:
    json.dump(data,f)

 plt.ioff()
 plt.show()
 

 #print(x)
