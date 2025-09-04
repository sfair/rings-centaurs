//-----------------------------------------------------------------------------------------//
//---- Simulação das partículas de um dos anéis de Saturno utilizando o integrador SEI ----//
//-----------------------------------------------------------------------------------------//
//---- José Eduardo e Rafael Sfair ----// 

//-------------------------------//
//---- Unidades da simulação ----//
//-------------------------------//
// Tempo -> [s]
// Velocidade -> [m/s]
// Distâncias (raio da partícula, tamanho da caixa, etc) -> [m]
// Massa -> [kg]
// Densidade Volumétrica (partículas) -> [kg/m^3]
// Densidade Superficial (anel) -> [kg/m^2]
// Constante Gravitacional -> [Nm^2/m^2]
// Frequência Epicíclica (Omega) -> [1/s]

//-----------------------------------------------//
//---- Importando as bibliotecas necessárias ----//
//-----------------------------------------------//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

//-------------------------------------------------------------//
//---- Definição de funções utilizadas durante a simulação ----//
//-------------------------------------------------------------//

// Adição do coeficiente de restituição uma vez que as velocidades deste exemplo dependem do mesmo
double coefficient_of_restitution_bridges(const struct reb_simulation* const r, double v)
{
    double eps = 0.32*pow(fabs(v)*100.,-0.234)  ;
    if (eps>1) eps=1                            ;
    if (eps<0) eps=0                            ;
    return eps                                  ;
}

// Declaração da função heartbeat
void heartbeat(struct reb_simulation* const r)  ;

//----------------------------------//
//---- Configurando a simulação ----//
//----------------------------------//
int main(int argc, char* argv[]) 
{
    struct reb_simulation* r = reb_simulation_create()                   ;
    
    
//---------------------------------------------------------------------------------------------------------------------------//
//---- Definição das constantes, do integrador e de métodos para resolução da gravidade, colisão e condições de contorno ----//
//---------------------------------------------------------------------------------------------------------------------------//
    r->opening_angle2     = .5                                           ;                   
    r->integrator         = REB_INTEGRATOR_SEI                           ;   
    r->boundary           = REB_BOUNDARY_SHEAR                           ;   
    r->gravity            = REB_GRAVITY_TREE                             ;     
    r->collision          = REB_COLLISION_TREE                           ;   
    r->collision_resolve  = reb_collision_resolve_hardsphere             ;  
    double OMEGA          = 0.00013143527                                ;        
    r->ri_sei.OMEGA       = OMEGA                                        ;                
    r->G                  = 6.67428e-11                                  ;         
    r->softening          = 0.1                                          ;                
    r->dt                 = 1e-3*2.*M_PI/OMEGA                           ;
    r->heartbeat          = heartbeat                                    ;        
    
//-----------------------------------------------//
//---- Definição das características do anel ----//
//-----------------------------------------------//
    double surfacedensity         = 400                                  ;        
    double particle_density       = 400                                  ;        
    double particle_radius_min    = 1                                    ;          
    double particle_radius_max    = 4                                    ;
    double particle_radius_slope  = -3                                   ;

//-----------------------------------------------//
//---- Definindo as características da caixa ----//
//-----------------------------------------------//
    double boxsize                = 100                                  ;
    reb_simulation_configure_box(r, boxsize, 2, 2, 1)                    ;
    r->N_ghost_x = 2                                                     ;
    r->N_ghost_y = 2                                                     ;
    r->N_ghost_z = 0                                                     ;
//É definido o tamanho da caixa (em x e y) e o número de root boxes (utilizadas para o paralelismo, se desejado)   
    
//------------------------------------//    
//---- Coeficiente de Restituição ----//
//------------------------------------//

//Atribui ao coeficiente de restituição o coeficiente de Bridges et al.
    r->coefficient_of_restitution = coefficient_of_restitution_bridges   ;
  
//Adição de uma velocidade de repulsão para que partículas com velocidade relativa nula não entrem uma na outra  
    r->minimum_collision_velocity = particle_radius_min*OMEGA*0.001      ; 

//-------------------------------------------//
//---- Adicionando as partículas do anel ----//
//-------------------------------------------//
    double total_mass = surfacedensity*r->boxsize.x*r->boxsize.y         ;
    double mass = 0                                                   ;
    while(mass<total_mass)
    {
        struct reb_particle pt;
        pt.x         = reb_random_uniform(r, -r->boxsize.x/2.,r->boxsize.x/2.)                                      ;
        pt.y         = reb_random_uniform(r, -r->boxsize.y/2.,r->boxsize.y/2.)                                      ;
        pt.z         = reb_random_normal(r, 1.)                                                                     ;
        pt.vx         = 0                                                                                           ;
        pt.vy         = -1.5*pt.x*OMEGA                                                                             ;
        pt.vz         = 0                                                                                           ;
        pt.ax         = 0                                                                                           ;
        pt.ay         = 0                                                                                           ;
        pt.az         = 0                                                                                           ;
        double radius     = reb_random_powerlaw(r, particle_radius_min,particle_radius_max,particle_radius_slope)   ;
        pt.r         = radius                                                                                       ;                                   

        double        particle_mass = particle_density*4./3.*M_PI*radius*radius*radius                              ;
        pt.m         = particle_mass                                                                                ;
        reb_simulation_add(r, pt)                                                                                   ;
        mass += particle_mass                                                                                       ;
    }
//As partículas são adicionadas de maneira aleatória dentro dos limites da caixa até que a massa total do anel seja atingida
//O raio das partículas obedece uma lei de potência aleatória dentro dos limites definidos
    
#ifdef OPENGL
    // Hack to artificially increase particle array.
    // This cannot be done once OpenGL is activated. 
    r->N_allocated *=8;
    r->particles = realloc(r->particles,sizeof(struct reb_particle)*r->N_allocated);
#endif // OPENGL

//Realizando a integração
    reb_simulation_integrate(r,100*(2*M_PI/OMEGA));
     // Cleanup
    reb_mpi_finalize(r);
    reb_simulation_free(r);
}

//-------------------------------------------------------------------------------//
//---- Determinação de funções para extrair quantidades físicas da simulação ----//
//-------------------------------------------------------------------------------//

//Construção da função responsável por medir a transparência do anel (profundidade óptica)
double mean_normal_geometric_optical_depth(const struct reb_simulation* const r)
{
    double area = 0.;
    for (int i=0;i<r->N;i++)
    {
        struct reb_particle p = r->particles[i];
        area += M_PI*p.r*p.r;
    }
    return area/(r->boxsize.x*r->boxsize.y);
}


double midplane_fillingfactor(const struct reb_simulation* const r){
    double area = 0.;
    for (int i=0;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        double R2 = p.r*p.r-p.z*p.z;
        if (R2>0.){
            area += M_PI*R2;
        }
    }
    return area/(r->boxsize.x*r->boxsize.y);
}

//Construção da estrutura responsável por calcular e atribuir a um array as velocidades de cada partícula
struct reb_vec3d velocity_dispersion(const struct reb_simulation* const r)
{
    struct reb_vec3d A = {.x=0, .y=0, .z=0}                                                     ;
    struct reb_vec3d Q = {.x=0, .y=0, .z=0}                                                     ;
    for (int i=0;i<r->N;i++)
    {
        struct reb_vec3d Aim1 = A                                                               ;
        struct reb_particle p = r->particles[i]                                                 ;
        A.x = A.x + (p.vx-A.x)/(double)(i+1)                                                    ;
        A.y = A.y + (p.vy+1.5*r->ri_sei.OMEGA*p.x-A.y)/(double)(i+1)                            ;
        A.z = A.z + (p.vz-A.z)/(double)(i+1)                                                    ;
        Q.x = Q.x + (p.vx-Aim1.x)*(p.vx-A.x)                                                    ;
        Q.y = Q.y + (p.vy+1.5*r->ri_sei.OMEGA*p.x-Aim1.y)*(p.vy+1.5*r->ri_sei.OMEGA*p.x-A.y)    ;
        Q.z = Q.z + (p.vz-Aim1.z)*(p.vz-A.z)                                                    ;
    }
    Q.x = sqrt(Q.x/(double)r->N)                                                                ;
    Q.y = sqrt(Q.y/(double)r->N)                                                                ;
    Q.z = sqrt(Q.z/(double)r->N)                                                                ;

    // Return velocity dispersion in xx, yy, zz
    return Q;
}
//O array A carrega as velocidades (não são retornadas) em relação ao shear e Q as velocidades em relação ao corpo central

double translational_viscosity(const struct reb_simulation* const r){
    double Wxy = 0.;
    for (int i=0;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        double vx = p.vx;
        double vy = p.vy+1.5*r->ri_sei.OMEGA*p.x;
        Wxy += vx*vy;
    }
    return 2./3.*Wxy/r->N/r->ri_sei.OMEGA;
}

double collisional_viscosity(const struct reb_simulation* const r){
    // This is a time average!
    // To reset, set r->collisions_plog equal to 0.
    double Mtotal = 0.;
    for (int i=0;i<r->N;i++){
        Mtotal += r->particles[i].m;
    }
    return 2./3./r->ri_sei.OMEGA/Mtotal/r->t* r->collisions_plog;
}


// Construção de uma função utilizada para salvar as posições de todas as partículas
// Como o objetivo é salvar as posições a cada período orbital os dados são salvos de forma dinâmica de acordo com a progressão da simulação
void save_positions(struct reb_simulation* const r, int period){
    char filename[50];
    sprintf(filename, "positions_T%d.aei", period);

    FILE* positions = fopen(filename, "w");

    for (int i = 0; i < r->N; i++) {
        fprintf(positions, "%e \t %e \t %e\n", r->particles[i].x, r->particles[i].y, r->particles[i].r);
    }
    
    fclose(positions);
}

//Função heartbeat chamada a cada passo de integração, salvando as velocidades e a profundidade ótica durante cada período orbital
//Um novo arquivo é criado assim que um período orbital é finalizado, salvando novamente todos os dados coletados durante o atual período
void heartbeat(struct reb_simulation* const r){
    //Definindo uma forma de contabilizar os períodos orbitais completos para que os dados sejam salvos da maneira correta
    static int period_count = 0;
    double orbital_period = 2.0 * M_PI / r->ri_sei.OMEGA;

    //Nomeando os arquivos de forma dinâmica e criando os mesmos
    char velocities_filename[50], optical_depth_filename[50];
    sprintf(velocities_filename, "velocities_T%d.txt", period_count + 1);
    sprintf(optical_depth_filename, "optical_depth_T%d.txt", period_count + 1);

    FILE* velocities = fopen(velocities_filename, "a");
    FILE* optical_depth_arq = fopen(optical_depth_filename, "a");

    //Coletando as informações de velocidade e profundidade ótica das funções que calculam essas grandezas
    struct reb_vec3d Q = velocity_dispersion(r);
    double optical_depth = mean_normal_geometric_optical_depth(r);

    fprintf(velocities, "%e\t%e\t%e\t%e\n", r->t, Q.x, Q.y, Q.z);
    fprintf(optical_depth_arq, "%e\t%e\n", r->t, optical_depth);

    fclose(velocities);
    fclose(optical_depth_arq);

    //Verificação de finalização do período orbital, caso tenha se completado as posições são salvas
    if (r->t >= (period_count + 1) * orbital_period) {
        period_count++;
    
        save_positions(r, period_count);
        printf("Arquivos referentes ao período %d gerados.\n", period_count);
    }
}
