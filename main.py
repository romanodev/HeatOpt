from js import Blob, document, URL
from pyscript import document
from heatopt import get_optimizer,create_stl,get_fig
from pyscript import display
import json
import numpy as np


structure = []


optimizer,grid = get_optimizer()

def run(event):

    #document.getElementById("loading").showModal()

    kxx = float(document.getElementById("input_kxx").value)
    kyy = float(document.getElementById("input_kyy").value)
    kxy = float(document.getElementById("input_kxy").value)


    #R   = document.getElementById("input_conic").value 
    if kxx > 0.0 and kxx < 1.0:
      if kyy > 0.0 and kyy < 1.0:
        if kxy > -0.1 and kxy < 0.1:
         
         kappa,fig,x,J = optimizer(kxx,kyy,kxy)
         structure[:] = x
         display(fig, target="chart",append=False)
         display(str(round(kappa[0],3)), target="output_kxx", append=False)
         display(str(round(kappa[1],3)), target="output_kyy", append=False)
         display(str(round(kappa[2],3)), target="output_kxy", append=False)
         #document.getElementById("loading").close()

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
 x = np.array(data['x'])    
 J = np.array(data['J'])    


 structure[:] = x

 fig = get_fig(x,grid,J)

 display(fig, target="chart")
 display(str(round(data['kappa'][0],3)), target="output_kxx", append=False)
 display(str(round(data['kappa'][1],3)), target="output_kyy", append=False)
 display(str(round(data['kappa'][2],3)), target="output_kxy", append=False)



