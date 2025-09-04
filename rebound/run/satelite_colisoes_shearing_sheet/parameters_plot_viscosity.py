import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from numpy.polynomial.polynomial import Polynomial
from matplotlib.ticker import FuncFormatter

#df = pd.read_csv(f"parameters_rank0_cf0.5.txt", sep='\t',header=None) 

# Dados fornecidos (já normalizados pelo fator comum)
factor = 0.00013143527 * 2.14**2
V_t = np.array([4.580406e-02,7.126021e-03,1.058478e-02,4.239076e-03,1.557025e-03,5.266498e-04,1.557757e-03]) / factor
#V_g = np.array([1.758906e-02,2.759955e-03,8.612841e-04,1.524198e-03,5.761390e-04,5.495190e-05]) / factor
V_col = np.array([9.700464e-03,3.600579e-03,1.443407e-03,7.272278e-04,3.675446e-04,2.001297e-04,9.002550e-04]) / factor
V_loc = np.array([1.851454e-02,7.654865e-04,8.280091e-03,1.987649e-03,6.133413e-04,2.715681e-04,6.575016e-04]) / factor

# Eixo x fictício
x = np.array([0.5,0.6,0.7,0.8,0.9,1.0,0.5])  # Substitua com seus valores reais

# Ajuste polinomial de grau 2 (pode mudar para 1 ou outro grau)

# Plot
plt.figure(figsize=(8, 6))

limy=np.linspace(0.05,25,10)
lim = np.full(10,0.9)
plt.plot(lim,limy,'--',alpha=0.6,color = 'purple',lw=0.8, label=r"$\tau$ = 0.8")

plt.plot(x, V_t, '--', label='total',color = 'black', lw=0.8)
plt.scatter(x, V_t, s=10,color = 'black')
'''
plt.plot(x, V_g, '--', label='gravitational',color = 'blue', lw=0.8)
plt.scatter(x, V_g, s=10,color = 'blue')

plt.plot(x, V_col, '--', label='collisional',color = 'grey', lw=0.8)
plt.scatter(x, V_col,s=10,color = 'grey')

plt.plot(x, V_loc, '--', label='local',color = 'red',lw=0.8)
plt.scatter(x, V_loc,s=10,color = 'red')
'''

#plt.ylim(0.05,25)
plt.xlabel(r'$r_h$')
plt.ylabel(r'$\nu$ $[\Omega r^2]$')
plt.xscale('log')
plt.yscale('log')


plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('viscosity_rh.png')

