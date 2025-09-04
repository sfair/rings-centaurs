/**
 * Shearing sheet with MPI
 * Giovana Ramon
 * 
 * This example simulates a small patch of Saturn's
 */

//--------Units---------//
// Teime -> [s]
// Velocity -> [m/s]
// Distances (Particle radius, boxsize, etc) -> [m]
// Mass -> [kg]
// Volumetric Density (particles) -> [kg/m^3]
// Superficial Density (ring) -> [kg/m^2]
// Gravitational Constant  -> [Nm^2/m^2]
// Epicic Frequency (Omega) -> [1/s]

//-------------------------------------------//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"


// ---------------------------------------------------------------------------------------------
// Restitution coefficient, constant value of 0.5

double coefficient_of_restitution_bridges(const struct reb_simulation* const r, double v)
{
    //double eps = 0.32*pow(fabs(v)*100.,-0.234);
    //if (eps>1) eps=1;
    //if (eps<0) eps=0;
    return 0.5;
}

//---------------------------------------------------------------------------------------------

void heartbeat(struct reb_simulation* const r)  ;

int main(int argc, char* argv[]) 
{
    struct reb_simulation* r = reb_simulation_create();
//---------------------------------------------------------------------------------------------

    //Start a time counter
    time_t inicio, fim;
    time(&inicio);

//---------------------------------------Rebound-----------------------------------------------

    r->opening_angle2     = .25;                   
    r->integrator         = REB_INTEGRATOR_SEI;   
    r->boundary           = REB_BOUNDARY_SHEAR;   
    r->gravity            = REB_GRAVITY_TREE;     
    r->collision          = REB_COLLISION_TREE;   
    r->collision_resolve  = reb_collision_resolve_hardsphere ;  
    double OMEGA          = 0.00013143527;   // epicic frequency to Saturn      
    r->ri_sei.OMEGA       = OMEGA;                                              
    r->G                  = 6.67428e-11;                               
    r->softening          = 0.1;                                   
    r->dt                 = 1e-3*2.*M_PI/OMEGA;                           
    r->heartbeat          = heartbeat;    
    r-> N_active          = 1;                                        
    
//-----------------Ring parameters---------------------------------------------------------------

    double surfacedensity         = 400;        
    double particle_density       = 400;   
    // Identical aprticles with an arbitrary radius                                    
    double particle_radius_min    = 1.46;                                              
    double particle_radius_max    = 1.46;                                 
    double particle_radius_slope  = 0;
    // Optical depth of the ring (this can be changed)
    double tau = 0.5;                                   
//---------------Sattelite parameters---------------------------
    double m_sat = 1000;
    double r_sat = 30;
//-------------------Box---------------------------//

    //double boxsize                = 100;  
    //Lambda Toomre definition                                
    double lambda_toomre =4.*M_PI*M_PI*surfacedensity/OMEGA/OMEGA*r->G;
    //Configure the boxsize in Toomre wavelength  
    double boxsize = 2.5 * lambda_toomre;
   
    if (argc>1){                      
        boxsize = atof(argv[1]);
    }

    // Setup 2x2 root boxes.
    // This allows you to use up to 4 MPI nodes.
    reb_simulation_configure_box(r, boxsize, 2, 2, 1);

    //---------------------MPI-----------------------//

    // Initialize MPI (this only works after reb_simulation_configure_box)
    reb_mpi_init(r);
    r->N_ghost_x = 2;
    r->N_ghost_y = 2;
    r->N_ghost_z = 0;                                                    

    r->coefficient_of_restitution = coefficient_of_restitution_bridges;

    // Repulsive velocity for particles with zero relative velocity do not enter each other 
    r->minimum_collision_velocity = particle_radius_min*OMEGA*0.001; 
      
    //------------------------Sattelites-----------//
      struct reb_particle sat = {0};
      sat.x = 0;//reb_random_uniform(r, -r->boxsize.x/4,r->boxsize.x/4);
      sat.y = 0;//reb_random_uniform(r, -r->boxsize.y/4,r->boxsize.y/4);

      sat.z = 0;
      sat.vx = 0;
      sat.vy = 0;
      sat.vz = 0;
      sat.ax = 0;
      sat.ay = 0;
      sat.az = 0;
      sat.m = m_sat;
      sat.r = r_sat;
      sat.hash = 1254863;
      //printf("P=%f\n", sat.period);
      reb_simulation_add(r,sat);

     // r->N_active = r->N;

#ifdef OPENGL
    // Hack to artificially increase particle array.
    // This cannot be done once OpenGL is activated. 
    r->N_allocated *=8;
    r->particles = realloc(r->particles,sizeof(struct reb_particle)*r->N_allocated);
#endif // OPENGL

//-------------------------------------------------------------------------------------------------------------

    //----------------- Start the integration-----------------------------//
    // Set the simulation time

    reb_simulation_integrate(r,10*2.*M_PI/OMEGA);

    // Cleanup
    reb_mpi_finalize(r);
    reb_simulation_free(r);

    // End the time counter
    for (long i = 0; i < 100000000; i++);
    time(&fim);
    double tempo_total = difftime(fim, inicio); 
    printf("Tempo decorrido: %.2f minutos\n", tempo_total / 60.0);
}

//--------------------------Functions to calculate physical quantities------------------//

//Optical Depth
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

//Filling factor
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

//Dispersion of velocities of the particles
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

// Local shear viscosity (acoording to Mondino & Salo 2019)
double translational_viscosity(const struct reb_simulation* const r){
    double Wxy = 0.;
    for (int i=0;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        double vx = p.vx;
        double vy = p.vy+1.5*r->ri_sei.OMEGA*p.x;
        Wxy += vx*vy;
    }
    return 2./3.*Wxy/r->N/r->ri_sei.OMEGA;   }

// Collisional viscosity (according to Mondino & Salo 2019)
double collisional_viscosity(const struct reb_simulation* const r){
    // This is a time average!
    // To reset, set r->collisions_plog equal to 0.
    double Mtotal = 0.;
    for (int i=0;i<r->N;i++){
        Mtotal += r->particles[i].m;
    }
    return 2./3./r->ri_sei.OMEGA/Mtotal/r->t* r->collisions_plog;
}

// Gravitational viscosity (according to Mondino & Salo 2019)
double gravitational_viscosity(const struct reb_simulation* const r){
    double G = r->G;
    double Omega = r->ri_sei.OMEGA;
    double torque = 0.;
    double Mtotal = 0.;

    // Sum of the torques between pairs (i<j)
    for (int i = 0; i < r->N; i++){
        Mtotal += r->particles[i].m;
    
        struct reb_particle pi = r->particles[i];
        for (int j = i+1; j < r->N; j++){
            struct reb_particle pj = r->particles[j];

            double dx = pj.x - pi.x;
            double dy = pj.y - pi.y;
            double dz = pj.z - pi.z;

            double r_ij2 = dx*dx + dy*dy + dz*dz;
            double r_ij = sqrt(r_ij2);
            double r3 = r_ij2 * r_ij ;

            double mprod = pi.m * pj.m;

            double torque_ij = mprod * dx * dy / r3;
            torque += torque_ij;
        }
    }

    double nu_grav = -2.0 * G / (3.0 * Omega * Mtotal) * torque;
    return nu_grav;
}



// Function to save the parameters in a file at each heartbeat
void heartbeat(struct reb_simulation* const r){
  if(reb_simulation_output_check(r,2.*M_PI/r->ri_sei.OMEGA)){   
    
      static int period_count = 0;
      //double orbital_period = 2.0 * M_PI / r->ri_sei.OMEGA;

     char parameters_filename[100];
     sprintf(parameters_filename, "parameters_rank%d_cf0.5.txt", r->mpi_id);

     FILE* parameters = fopen(parameters_filename, "a");

     struct reb_vec3d Q = velocity_dispersion(r);
     double optical_depth_val = mean_normal_geometric_optical_depth(r);
     double midplane_ff_val = midplane_fillingfactor(r);
     double transl_visc = translational_viscosity(r);
     double collis_visc = collisional_viscosity(r);
     double grav_visc = gravitational_viscosity(r);

     fprintf(parameters, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", 
            r->t, Q.x, Q.y, Q.z, 
            optical_depth_val, 
            midplane_ff_val,
            transl_visc,
            collis_visc,
            grav_visc,
            transl_visc+collis_visc+grav_visc);  

     fclose(parameters);

     struct reb_particle* sat_ptr = reb_simulation_particle_by_hash(r,1254863);
     FILE* sattelite = fopen("sattelite.txt","a");
     fprintf(sattelite, "%e\t%e\t%e\t%e\t%e\t%e\n", r->t,sat_ptr->x,sat_ptr->y,sat_ptr->z,sat_ptr->r,sat_ptr->m);
     fclose(sattelite);

     printf("Arquivos referentes ao período %d gerados.\n", period_count);
     period_count++;
     
      
}


}

