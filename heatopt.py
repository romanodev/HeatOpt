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

mpl.rcParams['toolbar'] = 'None'

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



def get_fig(x,grid):

  cmap = ListedColormap(["#cadbde","#404040"])

  if type(x) == list: x = np.array(x)  
  fig  = plt.figure(facecolor='#404040')
  ax   = fig.add_subplot(111)
  ax.imshow(1-np.pad(x.reshape((grid,grid)).T,grid,mode='wrap'),vmin=0,vmax=1,cmap=cmap,animated=True)
  ax.plot([grid,grid,2*grid,2*grid,grid],[grid,2*grid,2*grid,grid,grid],color='g',ls='--',lw=3)
  ax.axis('off')

  fig.tight_layout(pad=0,h_pad=0,w_pad=0)
  return fig


def compute_kappa_and_gradient(T,rho_dep,direction,N):

      dim = 2
      [kappa_ij,P_vec,i_mat,j_mat,ii,kappad_ij,kappad_ji] = rho_dep
      #Compute kappa---
      kappa     = np.array([np.sum(kappa_ij[ii[i]])*direction[i] - np.dot(T,P_vec[i]) for i in range(dim)])

      gradient = np.zeros(N)
      np.add.at(gradient,i_mat,kappad_ij*np.power(T[i_mat]-T[j_mat],2))
      for i in range(dim):
       tmp      = direction[i]**2 + 2*direction[i]*(T[j_mat[ii[i]]]-T[i_mat[ii[i]]])   
       np.add.at(gradient,i_mat[ii[i]],kappad_ij[ii[i]]*tmp)
       np.add.at(gradient,j_mat[ii[i]],kappad_ji[ii[i]]*tmp)

      return np.dot(kappa,direction),gradient


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



def get_aux(rho,i_mat,j_mat,ind_extremes):
    
     k0 = 1e-12;k1 = 1.0
     kappa_map       = k0 + rho*(k1-k0)
     kappa_ij     = 2*kappa_map[i_mat] * kappa_map[j_mat]/(kappa_map[i_mat] + kappa_map[j_mat])

     kappad_ij    = (k1-k0)*0.5*np.power(kappa_ij/kappa_map[i_mat],2)
     kappad_ji    = (k1-k0)*0.5*np.power(kappa_ij/kappa_map[j_mat],2)

     #assembly
     v1 = np.zeros_like(rho)
     np.add.at(v1,i_mat,kappa_ij)
     data = np.hstack((-kappa_ij,v1))
 
     N   = len(rho)
     dim = len(ind_extremes)

     i = np.hstack((i_mat,np.arange(N)))
     j = np.hstack((j_mat,np.arange(N)))

     P_vec   = np.zeros((dim,N))

     #compute vectorial perturbation
     for i in range(dim):
      ii  = ind_extremes[i,1]  
      np.add.at(P_vec,(i,i_mat[ii]), kappa_ij[ii])
      np.add.at(P_vec,(i,j_mat[ii]),-kappa_ij[ii])

     return data,i,j,kappa_ij,kappad_ij,kappad_ji,P_vec



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

    #For direct
    row = np.hstack((i_mat,np.arange(N))) 
    col = np.hstack((j_mat,np.arange(N))) 

    def func(rho):
     
     data,i,j,kappa_ij,kappad_ij,kappad_ji,P_vec = get_aux(rho,i_mat,j_mat,ind_extremes)

     aux = kappa_ij,i_mat,j_mat

     #Assemble matrix
     tmp = np.zeros_like(rho)
     np.add.at(tmp,i_mat,kappa_ij)
     d   = np.hstack((-kappa_ij,tmp))
     A   = scipy.sparse.csc_matrix((d, (row, col)), shape=(N, N))
     A[12,:] = 0;A[12,12] = 1
     lu = sp.linalg.splu(A)
     #---------------------

     #---------------------------------------
     kappa     = np.zeros(n_dir)
     jacobian  = np.zeros((n_dir,N))
     for n,direction in enumerate(directions):

      P = np.einsum('i,ic->c',direction,P_vec)   
     
      rho_dep = [kappa_ij,P_vec,i_mat,j_mat,ii,kappad_ij,kappad_ji]

      kappa_and_gradient =  partial(compute_kappa_and_gradient,rho_dep=rho_dep,N=N,direction=direction)

      #DIRECT-----------------------------------------
      P[12] = 0 
      T_mat  = lu.solve(np.array(P))
      kappa[n],jacobian[n] = kappa_and_gradient(T_mat)
      #-------------------------------------------------

     
     kappa    = np.matmul(Ainv,kappa)
     jacobian = np.einsum('ji,ik->jk',Ainv,jacobian)

     return kappa,jacobian

   
    return func


def get_optimizer():

 nb = 14
 mm = 30
 R  = 0.05
 grid    = 36  #Grid resolution
 N = grid**2
 betas = [2**(n+1) for n in range(nb)]
 maxiter = [mm]*nb
 maxiter.append(1)
 betas.append(1e24)

 f = fourier(grid=grid)

 mapping = get_mapping(grid,R=R)

 counts = [0]
 def optimize(*kd):
 
  kd = np.array(kd)  
  n_dir = len(kd)
  kappa = np.zeros_like(kd)
  def func(x,grad,beta):
   counts[0] +=1   
   x,gp =  mapping(x,beta)
   kappa[:],jacobian = f(x)
   g = np.linalg.norm(kappa-kd)
   jacobian = np.einsum('ik,kl->il',jacobian,gp)
   
   grad[:] = np.zeros(N)
   for i in range(n_dir):
    grad[:]  += 1/g*((kappa[i]-kd[i])*jacobian[i])

   return g

  opt = nlopt.opt(nlopt.LD_CCSAQ,N)

  x = np.random.rand(N)
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
  #print(np.sort(gradient_norm)[:30])
  x[gradient_norm < 1e-8] = 0
  #------------------

  fig = get_fig(x,grid)

  #print('total_count')
  #print(total_count)
  return kappa,fig,x

 return optimize,grid


if __name__ == '__main__':

 optimizer,grid = get_optimizer()

 kappa,fig,x = optimizer(0.1,0.2,-0.05)

 print(kappa)

 with open('structure.stl','w') as f:
  f.write(create_stl(x))

 plt.show()
