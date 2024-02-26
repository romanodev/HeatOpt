from js import Blob, document, URL
from pyscript import document
from heatopt import get_optimizer,create_stl,get_fig
from pyscript import display
import json

structure = []


optimizer,grid = get_optimizer()

def run(event):

    document.getElementById("loading").showModal()

    kxx = float(document.getElementById("input_kxx").value)
    kyy = float(document.getElementById("input_kyy").value)


    #R   = document.getElementById("input_conic").value 
    if kxx > 0.0 and kxx < 1.0:
      if kyy > 0.0 and kyy < 1.0:

        kappa,fig,x = optimizer(kxx,kyy)
        structure[:] = x
        display(fig, target="chart",append=False)
        display(str(round(kappa[0],3)), target="output_kxx", append=False)
        display(str(round(kappa[1],3)), target="output_kyy", append=False)
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


