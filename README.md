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

python heatopt.py
```





