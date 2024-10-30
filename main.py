from js import Blob, document, URL
from pyscript import document
from heatopt import get_optimizer,write_stl,get_fig
from pyscript import display
import json
from itertools import product
import numpy as np
from stl import mesh,Mode
from  io import BytesIO
import base64



structure = []


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


optimizer,grid = get_optimizer()

def run(event):

    document.getElementById("loading").showModal()

    kxx = float(document.getElementById("input_kxx").value)
    kyy = float(document.getElementById("input_kyy").value)
    kxy = float(document.getElementById("input_kxy").value)


    #R   = document.getElementById("input_conic").value 
    if kxx > 0.0 and kxx < 1.0:
      if kyy > 0.0 and kyy < 1.0:
        if kxy > -0.1 and kxy < 0.1:
         
         kappa,fig,x = optimizer(kxx,kyy,kxy)
         structure[:] = x
         display(fig, target="chart",append=False)
         display(str(round(kappa[0],3)), target="output_kxx", append=False)
         display(str(round(kappa[1],3)), target="output_kyy", append=False)
         display(str(round(kappa[2],3)), target="output_kxy", append=False)
         document.getElementById("loading").close()

def download_stl(event):
    """Create structure"""

    output = create_stl(structure)

    # content is the data to write to a file
    tag = document.createElement('a')
    blob = Blob.new([output], {type: "text/plain"})
    tag.href = URL.createObjectURL(blob)
    tag.download = 'structure.stl'
    tag.click()



if __name__ =='__main__':


 #Write example
 with open('./example.json', 'r') as fp:
     data = json.load(fp)


 structure[:] = data['x']

 fig = get_fig(data['x'],grid)

 display(fig, target="chart")
 display(str(round(data['kappa'][0],3)), target="output_kxx", append=False)
 display(str(round(data['kappa'][1],3)), target="output_kyy", append=False)
 display('0', target="output_kyy", append=False)



