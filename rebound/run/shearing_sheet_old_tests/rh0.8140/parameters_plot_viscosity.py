import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from numpy.polynomial.polynomial import Polynomial

df = pd.read_csv(f"parameters_rank0_cf0.5.txt", sep='\t',header=None) 

# Dados fornecidos (já normalizados pelo fator comum)
factor = 0.00013143527 * 2.14**2
V_t = np.array([3.613982e-04,4.031416e-04,1.026713e-03,2.076580e-03,4.745104e-03, 4.936365e-03, 1.075318e-02, 6.470153e-03, 9.285773e-03,1.281063e-02]) / factor
V_g = np.array([1.858497e-04,6.147289e-05,4.556897e-04,7.157111e-04,2.614604e-03, 2.298158e-03, 7.571851e-04, 7.864466e-04, 3.481772e-03,8.205951e-03]) / factor
V_col = np.array([4.038511e-05,9.939976e-05,1.954025e-04,3.612493e-04,6.526106e-04, 9.659412e-04, 1.331617e-03, 3.529189e-04, 1.617633e-03,1.100213e-03]) / factor
V_loc = np.array([1.351635e-04,2.422690e-04,3.756208e-04,9.996196e-04,1.477889e-03, 1.672266e-03, 8.664381e-03, 5.330788e-03, 3.546621e-03,3.504469e-03]) / factor

# Eixo x fictício
x = np.array([0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9,1.0])  # Substitua com seus valores reais

# Ajuste polinomial de grau 2 (pode mudar para 1 ou outro grau)

# Plot
plt.figure(figsize=(8, 6))

limy=np.linspace(0.05,25,10)
lim = np.full(10,0.8)
plt.plot(lim,limy,'--',alpha=0.6,color = 'purple',lw=0.8)

plt.plot(x, V_t, '--', label='total',color = 'black', lw=0.8)
plt.scatter(x, V_t, s=10,color = 'black')

plt.plot(x, V_g, '--', label='gravitational',color = 'blue', lw=0.8)
plt.scatter(x, V_g, s=10,color = 'blue')

plt.plot(x, V_col, '--', label='collisional',color = 'grey', lw=0.8)
plt.scatter(x, V_col,s=10,color = 'grey')

plt.plot(x, V_loc, '--', label='local',color = 'red',lw=0.8)
plt.scatter(x, V_loc,s=10,color = 'red')


plt.ylim(0.05,25)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\nu$ $[\Omega r^2]$')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('viscosity.png')

