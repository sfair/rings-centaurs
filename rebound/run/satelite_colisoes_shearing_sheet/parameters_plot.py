import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

for n in range(1, 11):  # loop sobre períodos (ajuste o range conforme necessário)
    all_pos = []  # lista para armazenar dados de todos os ranks
   
    for i in range(4):  # considerando 4 ranks (0 a 3)
        df = pd.read_csv(f"colisao_analise/parameters_T{n}_rank{i}_cf0.5.txt", sep='\t',header=None)
       
        df = df.astype(float)
        all_pos.append(df)  # adiciona os dados deste rank à lista
    
    # Concatenando os dados de todos os ranks em um único DataFrame
    pos = pd.concat(all_pos, ignore_index=True)

    
    tau = pos[4]
    FF = pos[5]
    V = pos[6] 
    vx = pos[1]
    vy = pos[2]
    vz = pos[3]

  
    plt.xlabel(r"$\tau$")
    plt.ylabel(r"vx")
    #plt.ylim(0,0.005)
    
    plt.scatter(tau,vx,s=0.3)

    plt.savefig(f'velocity.png')
    plt.close()

