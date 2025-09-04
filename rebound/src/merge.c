/**
 * @file    merge.c
 * @brief   Hard sphere collision resolution.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section     LICENSE

 */


#include <stdio.h>
#include <stdlib.h>
#include "particle.h"
#include "merge.h"
#include "rebound.h"
#include "boundary.h"
#include "tree.h"
#include "integrator_trace.h"
#ifdef MPI
#include "communication_mpi.h"
#endif // MPI
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

static void reb_tree_get_nearest_neighbour_in_cell(struct reb_simulation* const r, struct reb_vec6d gb, struct reb_vec6d gbunmod, int ri, double p1_r,  double second_largest_radius, struct reb_collision* collision_nearest, struct reb_treecell* c);
static void reb_tree_check_for_overlapping_trajectories_in_cell(struct reb_simulation* const r, struct reb_vec6d gb, struct reb_vec6d gbunmod, int ri, double p1_r, double p1_r_plus_dtv, struct reb_collision* collision_nearest, struct reb_treecell* c, double maxdrift);



int reb_collision_merge(struct reb_simulation* const r, struct reb_collision c){
    if (r->particles[c.p1].last_collision==r->t || r->particles[c.p2].last_collision==r->t) return 0;

    // Every collision will cause two callbacks (with p1/p2 interchanged).
    // Always remove particle with larger index and merge into lower index particle.
    // This will keep N_active meaningful even after mergers.
    int swap = 0;
    unsigned int i = c.p1;
    unsigned int j = c.p2;   //want j to be removed particle
    if (j<i){
        swap = 1;
        i = c.p2;
        j = c.p1;
    }

    struct reb_particle* pi = &(r->particles[i]);
    struct reb_particle* pj = &(r->particles[j]);

    double invmass = 1.0/(pi->m + pj->m);

    //Scale out energy from collision - initial energy
    double Ei=0, Ef=0;
    if(r->track_energy_offset){
        {
            double vx = pi->vx;
            double vy = pi->vy;
            double vz = pi->vz;
            // Calculate energy difference in inertial frame
            if (r->integrator == REB_INTEGRATOR_MERCURIUS && r->ri_mercurius.mode==1){
                vx += r->ri_mercurius.com_vel.x;
                vy += r->ri_mercurius.com_vel.y;
                vz += r->ri_mercurius.com_vel.z;
            }

            if (r->integrator == REB_INTEGRATOR_TRACE && r->ri_trace.mode==REB_TRACE_MODE_KEPLER){
                vx += r->ri_trace.com_vel.x;
                vy += r->ri_trace.com_vel.y;
                vz += r->ri_trace.com_vel.z;
            }

            Ei += 0.5*pi->m*(vx*vx + vy*vy + vz*vz);
        }
        {
            double vx = pj->vx;
            double vy = pj->vy;
            double vz = pj->vz;
            if (r->integrator == REB_INTEGRATOR_MERCURIUS && r->ri_mercurius.mode==1){
                vx += r->ri_mercurius.com_vel.x;
                vy += r->ri_mercurius.com_vel.y;
                vz += r->ri_mercurius.com_vel.z;
            }

            if (r->integrator == REB_INTEGRATOR_TRACE && r->ri_trace.mode==REB_TRACE_MODE_KEPLER){
                vx += r->ri_trace.com_vel.x;
                vy += r->ri_trace.com_vel.y;
                vz += r->ri_trace.com_vel.z;
            }

            Ei += 0.5*pj->m*(vx*vx + vy*vy + vz*vz);
        }
        const unsigned int N_active = ((r->N_active==-1)?r->N-r->N_var: (unsigned int)r->N_active);
        // No potential energy between test particles
        if (i<N_active || j<N_active){
            double x = pi->x - pj->x;
            double y = pi->y - pj->y;
            double z = pi->z - pj->z;
            double _r = sqrt(x*x + y*y + z*z);

            Ei += - r->G*pi->m*pj->m/_r;
        }
    }

    // Merge by conserving mass, volume and momentum
    pi->vx = (pi->vx*pi->m + pj->vx*pj->m)*invmass;
    pi->vy = (pi->vy*pi->m + pj->vy*pj->m)*invmass;
    pi->vz = (pi->vz*pi->m + pj->vz*pj->m)*invmass;
    pi->x  = (pi->x*pi->m + pj->x*pj->m)*invmass;
    pi->y  = (pi->y*pi->m + pj->y*pj->m)*invmass;
    pi->z  = (pi->z*pi->m + pj->z*pj->m)*invmass;
    pi->m  = pi->m + pj->m;
    pi->r  = cbrt(pi->r*pi->r*pi->r + pj->r*pj->r*pj->r);
    pi->last_collision = r->t;


    // Keeping track of energy offst
    if(r->track_energy_offset){
        {
            double vx = pi->vx;
            double vy = pi->vy;
            double vz = pi->vz;
            if (r->integrator == REB_INTEGRATOR_MERCURIUS && r->ri_mercurius.mode==1){
                vx += r->ri_mercurius.com_vel.x;
                vy += r->ri_mercurius.com_vel.y;
                vz += r->ri_mercurius.com_vel.z;
            }

            if (r->integrator == REB_INTEGRATOR_TRACE && r->ri_trace.mode==REB_TRACE_MODE_KEPLER){
                vx += r->ri_trace.com_vel.x;
                vy += r->ri_trace.com_vel.y;
                vz += r->ri_trace.com_vel.z;
            }

            Ef += 0.5*pi->m*(vx*vx + vy*vy + vz*vz);
        }
        r->energy_offset += Ei - Ef;
    }

    return swap?1:2; // Remove particle p2 from simulation
}