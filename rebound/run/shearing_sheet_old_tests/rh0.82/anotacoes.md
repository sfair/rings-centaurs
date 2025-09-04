## Anotações:
> **Shearing sheet 2 example:**\
  In this collision resolve method, particles are displaced if they overlap.\
`collision_resolve_hardsphere_pullaway`

>**Parâmetro rh:**
> $$r_h =  \frac{R_{Hill}}{2r_p}$$
> com $R_{Hill} = (\frac{2(M_1+M_2)}{3M_{planet}})^{\frac{1}{3}}a$\
> Here, $M_1$ and $M_2$ stand for the masses of the particles, and $M_{planet}$ stands for the mass of the central body. 


## Salvar em h5
- É possível salvar os dados do REBOUND em formato HDF5 de forma manual, utilizando 
 bibliotecas externas como h5py para python, exportando os dados relevantes.

[Caso do salvamento em Python com a biblioteca h5py](https://docs.h5py.org/en/stable/)

- Para C, é possivel importar dentro do makefile com:
`#include <hdf5.h>`

[Exemplo de uso em C](https://manual.nexusformat.org/examples/code_native.html
) 

[Biblioteca a implementação em C](https://support.hdfgroup.org/documentation/hdf5/latest/_u_g.html)


### To do (Reunião)

- reduzir o salvamento pra cada período e não passo de integração 
> if para salvamento\
> `reb_simulation_output_check`

- Incrementar a viscosidade local e gravitacional no código

- Incrementar o rh:
> rh was defined as the ratio of the Hill radius for a pair
of identical particles over the sum of their physical radii

- Colocar os aquivos no github
- Checar a viscosidade e os parâmetros rh

> Importante!
![alt text](image.png)