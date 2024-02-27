# HeatOpt

HeatOpt is a "computational paper" showcasing a heat transport topology optimization solver accompanied by technical details. 

The main goal is to optimize a material with a prescribed effective thermal conductivity tensor of the form

$$\kappa = \begin{pmatrix}
\kappa_{xx} & 
\kappa_{xy} \\
\kappa_{yx} & 
\kappa_{xy}
\end{pmatrix}.$$


The engine is written in Python and is based on [PyScript](https://pyscript.net/)

Author: Giuseppe Romano (romanog@mit.edu)

Computational Paper: [Link](https://romanodev.github.io/HeatOpt/)

To run HeatOpt as a stand-alone program, run the following code

```bash
python -m venv env
source env/bin/activate
pip install -r requirements.txt
```

Optimizing a material with tensor

$$\kappa = \begin{pmatrix}
0.1 & 
-0.04 \\
-0.04 & 
0.2
\end{pmatrix}. $$

can be achieved with the following code

```python

from heatopt import get_optimizer,write_stl
import matplotlib.pylab as plt


optimizer,grid = get_optimizer()

kappa,fig,x = optimizer(0.1,0.2,-0.04)

print(kappa)

with open('structure.stl','w') as f:
  f.write(write_stl(x))

plt.savefig('structure.png')
plt.show()
```

![Alt text](structure.png)




