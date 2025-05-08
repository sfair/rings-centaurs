import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

V_t= [4.745104e-03/ (0.00013143527*2.14**2),4.936365e-03/ (0.00013143527*2.14**2),1.075318e-02/ (0.00013143527*2.14**2),6.470153e-03/ (0.00013143527*2.14**2)]
V_g = [2.614604e-03/ (0.00013143527*2.14**2),2.298158e-03/ (0.00013143527*2.14**2),7.571851e-04/ (0.00013143527*2.14**2),7.864466e-04/ (0.00013143527*2.14**2)]
V_col = [6.526106e-04/ (0.00013143527*2.14**2),9.659412e-04/ (0.00013143527*2.14**2),1.331617e-03/ (0.00013143527*2.14**2),3.529189e-04/ (0.00013143527*2.14**2)]
V_loc = [1.477889e-03/ (0.00013143527*2.14**2),1.672266e-03/ (0.00013143527*2.14**2),8.664381e-03/ (0.00013143527*2.14**2),5.330788e-03/ (0.00013143527*2.14**2)]
    #tau = df[4]
tau = [0.5,0.6,0.7,0.8]
FF = df[5]
V = df[6]
vx = df[1]
vy = df[2]
vz = df[3]
V_total = df[9]/ (0.00013143527*2.14**2)

  
#plt.xlabel(r"$\tau$")
plt.ylabel(r"$FF$")
plt.ylabel(r"$\nu [\Omega r^2]$")
    #plt.ylim(0,0.005)
#tau = np.full(len(V_t),0.5)
plt.scatter(tau,V_t, label="Total viscosity")
plt.scatter(tau,V_g, label="Gravitational viscosity")
plt.scatter(tau,V_col,  label="Collisional viscosity")
plt.scatter(tau,V_loc,label="Local viscosity")
plt.xscale('log')
plt.yscale('log')
plt.legend()


plt.savefig(f'viscosity.png')
plt.close()

