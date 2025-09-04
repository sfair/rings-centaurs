/**
 * @file    hardsphere.c
 * @brief   Hard sphere collision resolution (colisão de esferas rígidas).
 * 
 * Resolve colisões detectadas, assumindo que as partículas são esferas rígidas.
 * Faz rotações para alinhar a linha de centros com o eixo x, resolve o choque
 * como um problema 1D com restituição, e depois desfaz a rotação.
 */

#include <stdio.h>
#include <stdlib.h>
#include "particle.h"
#include "hardsphere.h"
#include "rebound.h"
#include "boundary.h"
#include "tree.h"
#include "integrator_trace.h"
#ifdef MPI
#include "communication_mpi.h"
#endif // MPI

// Macros auxiliares
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< mínimo
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< máximo


int reb_collision_resolve_hard(struct reb_simulation* const r, struct reb_collision c){
    // Ponteiro para o array de partículas
    struct reb_particle* const particles = r->particles;

    // Cópias locais de p1 e p2 (não serão alteradas diretamente)
    struct reb_particle p1 = particles[c.p1];
    struct reb_particle p2 = particles[c.p2];

    // Vetor "ghost box" (correção para condições periódicas)
    struct reb_vec6d gb = c.gb;

    // Vetor posição relativa entre os centros (já corrigido por periodicidade)
    double x21  = p1.x + gb.x  - p2.x; 
    double y21  = p1.y + gb.y  - p2.y; 
    double z21  = p1.z + gb.z  - p2.z; 

    // Distância de contato (soma dos raios)
    double rp   = p1.r+p2.r;

    // Guardar vy da partícula "externa" (para log depois)
    double oldvyouter;
    if (x21>0){
        oldvyouter = p1.vy;
    }else{
        oldvyouter = p2.vy;
    }

    // Teste de não contato: se distância > soma dos raios, não há colisão
    if (rp*rp < x21*x21 + y21*y21 + z21*z21) return 0;

    // Velocidade relativa (corrigida por periodicidade)
    double vx21 = p1.vx + gb.vx - p2.vx; 
    double vy21 = p1.vy + gb.vy - p2.vy; 
    double vz21 = p1.vz + gb.vz - p2.vz; 

    // Teste de aproximação: produto escalar v·r > 0 significa que estão se afastando
    if (vx21*x21 + vy21*y21 + vz21*z21 >0) return 0;

    // ======== ROTACIONAR SISTEMA ========

    // Rotação em torno de x para colocar vetor no plano xy
    double theta = atan2(z21,y21);
    double stheta = sin(theta);
    double ctheta = cos(theta);
    double vy21n = ctheta * vy21 + stheta * vz21;    
    double y21n  = ctheta * y21  + stheta * z21;    

    // Rotação em torno de z para alinhar vetor com eixo +x
    double phi = atan2(y21n,x21);
    double cphi = cos(phi);
    double sphi = sin(phi);
    double vx21nn = cphi * vx21  + sphi * vy21n;        

    // ======== RESOLVER COLISÃO 1D NA DIREÇÃO NORMAL ========

    // Coeficiente de restituição (default = 1, elástico)
    double eps= 1;
    if (r->coefficient_of_restitution){
        eps = r->coefficient_of_restitution(r, vx21nn);
    }

    // Variação da velocidade relativa normal após colisão
    double dvx2 = -(1.0+eps)*vx21nn;

    // Impulso mínimo para evitar que partículas "grudem"
    double minr = (p1.r>p2.r)?p2.r:p1.r;
    double maxr = (p1.r<p2.r)?p2.r:p1.r;
    double mindv= minr*r->minimum_collision_velocity;
    double _r = sqrt(x21*x21 + y21*y21 + z21*z21);
    mindv *= 1.-(_r - maxr)/minr;
    if (mindv>maxr*r->minimum_collision_velocity) mindv = maxr*r->minimum_collision_velocity;
    if (dvx2<mindv) dvx2 = mindv;

    // ======== DESFAZER ROTAÇÕES (voltar ao sistema original) ========
    double dvx2n  = cphi * dvx2;        
    double dvy2n  = sphi * dvx2;        
    double dvy2nn = ctheta * dvy2n;    
    double dvz2nn = stheta * dvy2n;    

    // ======== APLICAR O IMPULSO NAS DUAS PARTÍCULAS ========
    // Repartição proporcional às massas, garantindo conservação de momento

    const double p2pf = p1.m/(p1.m+p2.m);
    particles[c.p2].vx -= p2pf*dvx2n;
    particles[c.p2].vy -= p2pf*dvy2nn;
    particles[c.p2].vz -= p2pf*dvz2nn;
    particles[c.p2].last_collision = r->t;

    const double p1pf = p2.m/(p1.m+p2.m);
    particles[c.p1].vx += p1pf*dvx2n; 
    particles[c.p1].vy += p1pf*dvy2nn; 
    particles[c.p1].vz += p1pf*dvz2nn; 
    particles[c.p1].last_collision = r->t;

    // ======== LOG DE COLISÕES (diagnóstico) ========
    if (x21>0){
        r->collisions_plog += -fabs(x21)*(oldvyouter-particles[c.p1].vy) * p1.m;
        r->collisions_log_n ++;
    }else{
        r->collisions_plog += -fabs(x21)*(oldvyouter-particles[c.p2].vy) * p2.m;
        r->collisions_log_n ++;
    }

    //if (c.p1 == 0){
      //  return 2;
    //}else if (c.p2 == 0){
      //  return 1;
    //}
        return 0; // colisão resolvida, nenhuma fusão
}
