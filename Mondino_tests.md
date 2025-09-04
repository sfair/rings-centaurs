## teste 0
- Padrão do shearing sheet com o coeficiente de restituição dependente da velocidade e o 
ângulo de abertura igual a 0.5
- Número de partículas = 5949
- Após esse teste mudei o salvamento dos arquivos para salvar todos os parâmetros calculados em um unico arquivo
e tornei o **coeficiente de restituição constante** 

## teste1 

- Ângulo de abertura modificado para 0.25, correspondendo com o artigo de Mondino & Salo, 2022
- Coeficiente de restituição igual a 0.5
- N = 6082

## teste 2

- Adicionei o comprimento de Toomre, sendo que o tamanho da caixa é igual a 5 $\lambda_{crit}\times\lambda_{crit}$

-  $\lambda_{crit}$ = 61.009899 m
- Rodando a simulação por 10 períodos orbitais 
- N = 13920

> **Anotação:**  
> Para ver informações dos cores:  
> `cat /proc/cpuinfo | less`



## 05/05 Testes: 

### teste 1:
- Com a nova forma de salvamento realizei um teste com 15 períodos em um único nodo, levou 44 minutos gerou 13837 partículas

- A profundidade óptica não sofreu variação a cada intervalo de período

> Profundiade óptica média:
> $$ \tau = \frac{N \pi r^2}{LxLy}$$

- As equaçãoes da viscosidade colisional e transversal batem com a viscosidade colisional e local presentes no artigo

- Viscosidade translacional: representa o transporte de momento angular radial mediado pelas flutuações de velocidade.
- $Wxy$: fluxo de momento angular específico por unidade de área, transportado pelo movimento aleatório das partículas. 

 ⚠️ **Tarefa Futura:** Implementar a viscosidade gravitacional


### Tentativa de rodar sem gravidade mútua para reproduzir os primeiros gráficos do artigo
- $\tau = 0.5$
- 10 períodos
- ⚠️ Não continuei com essas modificações 

## Ajuste o codigo para a variação do $\tau$ a cada simulação e numero de particulas 

- Primeiro teste:

         Parâmetros

         boxsize = 2 $\lambda_{crit}$
         rh = 0.81
         rp = 1.18m
         Densidade da partícula = 400 kg/m^3
         Densidade superficial do anel = 400kg/m^2

    - **$\tau$ = 0.5** 
         -  10 periodos 
        - N = 6807
    
        - Duração: 20 min

    - **$\tau$ =1.0**
        - Number of particules = 13614
        - Tempo decorrido: 60.30 minutos

    - **$\tau$ = 0.6** 
        - Number of particules = 8168
        - Tempo decorrido: 16.47 minutos



## 07/05

### Caso de Saturno
> Quem define o rh nessa simulação vai
> ser o $\Omega = \sqrt{\frac{GM}{a^3}}$,
> uma vez que alteramos o semi eixo, mudamos o omega da simulação e assim o rh. Isso acontece porque rh depende do Raio de hill e do raio da partícula, se mantemos constantes a densidade da partícula, a densidade superficial e o raio da partícula, a única maneira de alterar o rh se torna o semi-eixo.

        Parâmetros:
        - rh = 0.81
        - r = 1.46 m
        - $\lambda_{crit}$ = 61.009899 m
        - $\rho$ = 400 kg/m³
        - Densidade superficial = 400 kg/m²
        - R_hill = 2.26m

------------------------------------------
-  $\tau$ = 0.5
    - Number of particules = 6947
    - 12.22 minutos
    - A ordem da viscosidade bateu! 

> ⚠️ Para resolver o problema da viscosidade bater com os resultados do artgio (no artigo $\nu [\Omega r^2]$ | na simulação $\nu [m^2/s]$) dividi o resultado na importação do plot por $\Omega$ e $r^2$, sendo que 
a unidade de saída do código é m^2/s


-------------------------------------------

-  $\tau$  = 0.6
    - Number of particules = 8337
    > Rodei no computador da feg em um nodo e na workstation em 4 nodos
    - **Resultado da workstation**: Tempo decorrido = 4.13 minutos
    - **Resultado para o computador da feg**: Tempo decorrido = 16.37 minutos
    - com MPI no computador da feg: Tempo decorrido: 41.58 minuto ❓❓❓

---------------------------------------------
-  $\tau$ = 0.7
    - Number of particules = 9724
    - Na Moria com 4 nodos: Tempo decorrido: 5.38 minutos

>  Não foi o resultado que eu usei por conta do problema com o MPI

---------------------------------------------

- $\tau$ = 0.8
    - Number of particles = 11116
    - Tempo decorrido: 32.25 minutos

<div style="border:1px solid black; padding:10px; background-color:#f0f8ff; color:darkblue;">
<span style="font-weight:bold;">Problema:</span> Quando rodo em MPI em 4 nodos diferentes estou tendo um problema com a viscosidade colisional que fica negativa em alguns casos, o que não ocorre quando não faço a separação.
<span style="font-weight:bold;">Função caixa preta na viscosidade colisional: collision_plog</span>
</div>

<div style="border:1px solid black; padding:10px; background-color:#f0f8ff; color:darkblue;">
<span style="font-weight:bold;">Conclusão:</span>  Realmente tenho um problema para calcular globalmente a viscosidade e talvez outros parâmetros que dependam 
da interação entre bordas na paralelização, o último ishue no git do hanno é sobre o mesmo problema 

### Entender o problema já resolvido e verificar a aplicabilidade para o meu caso. 
</div>


[Rebound issue](https://github.com/hannorein/rebound/blob/7d7594760072cb60d6d42a2b2f8c6487d2ee077f/examples/mpi_unittests/problem.c)

> Continuando as simulações para um único nodo

- $\tau$ = 0.9
    - Number of particules = 12506
    - Tempo decorrido: 41.17 minutos

- $\tau$ = 1.0
    - Number of particules = 13895
    - Tempo decorrido: 50 minutos


## 08/05

> Vou fazer as simulações com tau entre 0.1 e 0.4 para completar o gráfico para o rh 0.81

- $\tau = 0.1$
    - Number of particules = 1389
    - Tempo decorrido: 1.50 minutos

- $\tau = 0.2$
    - Number of particules = 2779
    - Tempo decorrido: 3.58 minutos

- $\tau = 0.3$
    - Number of particules = 4168
    - Tempo decorrido: 6.50 minutos

- $\tau = 0.4$
    - Number of particules = 5558
    - Tempo decorrido: 8.12 minutos
------------------------------------------------

> Verificando o semi eixo para rh = 0.82\
> **Resultado**\
> a = 1.30923403335e+8\
> Omega = 0.000130007


> **Para o omega padrão do shearing sheet** \
> rh = 0.8140\
> a = 1.2997e8  m

----------------------------------------------

### Novos casos com rh =0.82
> Parâmetros:\
> a = 1.309e8m\
> omega = 0.000130007 m²/s\
> r = 1.46m\
> densidades = 400 kg/m² e kg/m³\
> Toomre wavelength: 61.009899m

- $\tau=0.1$
 
    - Number of particules = 1389
    - Tempo decorrido: 1.52 minutos

- $\tau = 0.2$
    - Number of particules = 2779
    - Tempo decorrido: 3.53 minutos


- $\tau = 0.3$
    - Number of particules = 4168
    - Tempo decorrido: 5.63 minutos

- $\tau = 0.4$
    - Number of particules = 5558
    - Tempo decorrido: 8.15 minutos

- $\tau =0.5$
    - Number of particules = 6947
    - Tempo decorrido: 11.02 minutos

- $\tau =0.6$
    - Number of particules = 8337
    - Tempo decorrido: 16 minutos

- $\tau = 0.7$
    -Number of particules = 9727
    - 17.78 minutos

- $\tau = 0.8$
    -Number of particules = 11116
    - 15 minutos
    - Meu computador desligou e foi só até 7 períodos

- $\tau = 0.9$
    - Number of particules = 12506
    - 39.73 minutos
    

> <span style="color:red;">Importante</span>\
> Vou ter que fazer novos testes mudando agora diretamente o valor do raio da particula no código de acordo com o $\tau$ e oraio da partícula, segundo a equação presente em Mondino & Salo 2025
> Já implementei isso nó código disponível na workstation Moria


