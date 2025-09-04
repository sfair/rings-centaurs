**não esquecer**: para rodar o make com MPI precisa que o open mpi esteja instalado no computador

## Shearing sheet com satélite 

- Inicialmente estava adicionando o satélite como uma partícula além das que já eram definidas pela profundidade óptica, entretanto o satélite se deslocava do centro devido a interação gravitacional mútua. 

### Teste simulando somente o satélite: 

- Removi a caixa de partículas do shearing sheet e deixei N_active = 1, sendo a massa do satélite igual a 1000kg e o raio igual a 30m em 10 períodos; 

- O satélite permaneceu no centro como esperado; 

### Teste com as partículas sem massa

- Adicionei novamente as partículas com raio 1.46 m, a massa do satélite 100kg e raio 50m , com a caixa de 2.5 $\lambda _c$, N_active comentado , os satélites adicionados primeiro no centro da caixa, ou seja, com as coordendas zeradas.
- Para salvar seu hash criei uma variável sat_index -> N 
- As partículas foram para $\tau = 0.5$ e  não tinham massa, m = 0kg. O período foi igual a 4. 

- O satélite permanceu no centro da caixa sendo que conforme os períodos evoluiram as particulas se dispersram ao redor do corpo.

### Teste com cavidade

- Para esse caso mantive praticamente os mesmos parâmetros do anterior, entretanto, com $\tau =0.8$, a massa ainda igual a zero para os satélites. 

- Para criar a cavidade foi feito: 

```c
double d2 = pt.x*pt.x + pt.y*pt.y + pt.z*pt.z 
if (d2<(1.1*r_sat)*(1.1*r_sat)){
    continue;
}
```

- Isso cria uma cavidade ao redor do satélite com tamanho equivalente ao raio igula 1.1 o raio de satélite, adicionando as partículas ao redor desse espaço sem exclui-las. 

- O satélite permaneceu no centro, as partículas ficaram ao redor da cavidade e conforme os períodos evoluiram até 10, foi se criando um espaço vazio em cima e embaixo do satélite, sendo que no último período quase não existem partículas.


## Remoção das partículas que encontram o satélite

### Criação de rotina de colisão personalizada no src
- Podemos criar nossa rotina personalizada de colisão adicionando no início do código uma função declarada `int reb_collision_resolve_status(struct reb_simulation* const r, struct reb_collision c);` e então chama-la no `r->collision_resolve  = reb_collision_resolve_status;  `

- No caso estudado estamos utilizando como busca de colisões `r->collision          = REB_COLLISION_TREE; `

    >This method uses an oct-tree to check for overlapping particles at the end of the timestep. When a large number of particles is used, this method scales as , rather than for the direct search. Note that you need to initialize the simulation box whenever you want to use the tree. Below is an example on how to enable the tree based collision search.C

- Para modificar funções no src podemos criar um arquivo .c dentro da pasta src e então criar um arquivo.h com o mesmo nome contendo:

    ```c
    /**
     * @file 	hardsphere.h
     * @brief 	Tools for creating distributions.
     * @author 	Hanno Rein <hanno@hanno-rein.de>
     * 
     * @section 	LICENSE
     * Copyright (c) 2013 Hanno Rein
     */
    #ifndef HARDSPHERE_H
    #define HARDSPHERE_H

    #endif
    ``` 
- E então incluir esse arquivo no `rebound.c` com `#include "hardsphere.h"` e adicionar ao ŕebound.h a função presente no arquivo.c DLLEXPORT `int reb_collision_resolve_hard(struct reb_simulation* const r, struct reb_collision c);`

- No Makefile do src adicionar o nome do arquivo n seção:

    ```c
    SOURCES=rebound.c tree.c particle.c gravity.c integrator.c integrator_whfast.c integrator_whfast512.c integrator_saba.c integrator_ias15.c integrator_sei.c integrator_bs.c integrator_leapfrog.c integrator_mercurius.c integrator_trace.c integrator_eos.c boundary.c input.c binarydiff.c output.c collision.c communication_mpi.c display.c tools.c rotations.c derivatives.c simulationarchive.c glad.c integrator_janus.c transformations.c fmemopen.c server.c hardsphere.c merge.c
    ```

### Aplicação no caso do shearing sheet com satélite

- Adicionei ao final do problem.c a função:

```c
int reb_collision_resolve_status(struct reb_simulation* const r, struct reb_collision c){
    int result ;

    struct reb_particle* p1 = &(r->particles[c.p1]);
    struct reb_particle* p2 = &(r->particles[c.p2]);

    if (p1->hash == SAT_HASH){
        result = 2;   // remove p2
    }else if (p2->hash == SAT_HASH){
        result = 1;   // remove p1
    }else{
        result = reb_collision_resolve_hard(r, c);
    }

    if (result != 0){
        col_sat++; 
        FILE* of_final = fopen("collision.csv", "a");
        if (of_final){
            fprintf(of_final, "%d\t%d\t%d\t%e\n", col_sat, c.p1, c.p2, r->t);
            fclose(of_final);
        }
    }

    col++;  
    return result;  
}
```
- Onde SAT_HASH é o hash do satélite definido globalmente antes do int mais, col_sat é o número de colisões com o satélite e col é o número total de colisões, definidos como contadores.

- O rebound tem como definição para o retorno da função de colisão:
    - 0: nada acontece
    - 1: remove p1
    - 2: remove p2
    - 3: remove ambos

- Assim, na função personalizada, verificamos se uma das partículas é o satélite e então removemos a outra partícula. Caso contrário, a colisão é resolvida como uma colisão elástica (hard sphere).

- As partículas estão sendo acessadas através do índice presente na estrutura de colisão `struct reb_collision c`, que contém os índices das partículas que colidiram, `c.p1` e `c.p2`.

> Parâmetros utilizados:

- Teste com condições de Saturno
   - OMEGA  = 0.00013143527 [1/s];
   - surfacedensity = 400 [kg/m^2];
   - particle_density  = 400 [kg/m^3];
   - particle_radius   = 1.46 [m];
   - tau = 0.5;      
   - boxsize = 2.5 * lambda_toomre
   - Number of particules = 6947
    - Satélite:
    - m_sat = 30159289.474 Kg; #correspondente a densidade de 0.9 g/cm³ e 20 m de raio
    - r_sat = 20 m;

> Partículas com massa zero
- As partículas foram adicionadas com massa zero para que a interação gravitacional mútua entre elas e o satélite fosse desprezada, evitando que o satélite se deslocasse do centro da caixa.

    Número total de colisões: 8542

    Número total de colisões com o satélite: 381

> Partículas com massa
- As partículas foram adicionadas com massa correspondente a densidade de 400 kg/m³ e raio de 1.46 m, o satélite se desloca um pouco em y limpando a área ao redor, mas permanece no centro em x.

    Número total de colisões: 5644083

    Número total de colisões com o satélite: 3520

### Teste aumenando o tamanho da caixa
- Mantive os mesmos parâmetros do teste anterior, mas aumentei o tamanho da caixa para 1000km em y e diminiu $\tau$ para 0.4

   - Boxsize_y = 1067.67 m
   - Boxsize_x = 305.05 m
   - Number of particules = 19454

- Diminui o passo de salvamento para 1/5 do período orbital


> A fazer: 
- Mudar a identificação da partícula a ser removida para verificar a quem tem maior massa, assim sempre vai identificar os atélite sem precisar forçar o hash. 

- Testar o MPI: mudar o makefile e adicionar as funções no problem.c
- Primeiro teste apresentou problemas com a remoção das partícualas em paralelo.
