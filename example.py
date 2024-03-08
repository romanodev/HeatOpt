from heatopt import get_optimizer,write_stl
import matplotlib.pylab as plt


optimizer,grid = get_optimizer()

kappa,fig,x = optimizer(0.4,0.3,-0.04)

print(kappa)

quit()
with open('structure.stl','w') as f:
  f.write(write_stl(x))

plt.savefig('structure.png')
plt.show()
