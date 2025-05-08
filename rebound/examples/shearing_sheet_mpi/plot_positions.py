import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# Number of periods
N = 11
# Number of ranks in the MPI
rank = 1

for n in range(1, N):  
    all_pos = []  # to save datas in all ranks
    
    for i in range(rank):  
        df = pd.read_csv(f"positions_T{n}_rank{i}_cf0.5.aei", delim_whitespace=True, header=None, usecols=[0,1,2])
        df = df.astype(float)
        all_pos.append(df)  
        
    
    pos = pd.concat(all_pos, ignore_index=True)

    X = pos[0].tolist()  
    Y = pos[1].tolist()  
    R = pos[2].tolist()  

    fig = plt.figure(figsize=(8,8))
    ax = plt.subplot(111, aspect='equal')
    ax.set_ylabel(r"y [m]")
    ax.set_xlabel(r"x [m]")
    ax.set_ylim(-2.5*61.009, 2.5*61.009)
    ax.set_xlim(-2.5*61.009, 2.5*61.009)    

    for x, y, r in zip(X, Y, R):   
        circ = patches.Circle((x, y), r, facecolor='darkgray', edgecolor='black')
        ax.add_patch(circ)      

    plt.savefig(f'particles_T{n}_N{N}_cf0.5.png')
    plt.close(fig)

