# HeatOpt

HeatOpt is a "computational paper" showcasing a heat transport topology optimization solver accompanied by technical details. 

The main goal is to optimize a material with a prescribed effective thermal conductivity tensor of the form

$$\kappa = \begin{pmatrix}
\kappa_{xx} & 
\kappa_{xy} \\
\kappa_{yx} & 
\kappa_{xy}
\end{pmatrix}.$$

The tensor is reconstructed by running three separate simulations, one for each direction of the applied temperature gradient. The thermal conductivity along a given direction $$\mathbf{\hat{n}}$$ is given by $$\kappa(\mathbf{\hat{n}}) = \mathbf{\hat{n}}^T\kappa \mathbf{\hat{n}} = n_x^2 \kappa_{xx}+ n_y^2 \kappa_{yy} + 2n_xn_y\kappa_{xy}$$.

We solve the heat conduction for the directions $$\mathbf{\hat{n}}_0 = \mathbf{\hat{x}}$$ , $$\mathbf{\hat{n}}_1 = \mathbf{\hat{y}}$$ , and $$\mathbf{\hat{n}}_2 = \sqrt{2}/2  \mathbf{\hat{x}} + \sqrt{2}/2  \mathbf{\hat{y}}  $$, leading to the linear system

```math
   \begin{bmatrix}
            \kappa_{xx}  \\
            \kappa_{yy} \\
            \kappa_{xy}  \\
          \end{bmatrix}   = \begin{bmatrix}
            1 & 0 & 0 \\
            0 & 1 & 0 \\
            -1/2 & -1/2 & 1 \\
          \end{bmatrix}
\begin{bmatrix}
            \kappa_0  \\
            \kappa_1 \\
            \kappa_2  \\
          \end{bmatrix} 
```



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
#Always check if the final kappa is close to the prescribed one
print(kappa)

with open('structure.stl','w') as f:
  f.write(write_stl(x))

plt.savefig('structure.png')
plt.show()
```
Note that the thermal conductivity for the base material is 1 W/m/K. 



![Alt text](structure.png)




