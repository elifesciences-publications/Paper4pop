#ifndef __GLOBALVARS__ 
#define __GLOBALVARS__ 

// Simulation time constants
#define DT .01
#define duration 10.0E3 // 10.E3  //

#define Tst 2.0E3 // 10.E3 // 
#define Tw 1.0E3 // 1.00E3 // 10.E3 // 
#define Tl 2.0E3 // 2.0E3 // 30.0E3 // 

// 10971.5
// #define nbPref 10000
#define nbPref 10000
#define IF_Nk 0

#define IF_EXP_DECAY 1 
#define IF_EXP_RISE_DECAY 0 

#define IF_BENCHMARK 0 
#define IF_INTERPOLATION 1
#define IF_EULER 0 
#define IF_RK2 1 

#define IF_EtoSOM_O1 0

//////////////////////////
// LIF parameters 
//////////////////////////
#define Vr -70. // Membrane Resting Potential 
#define Vth -50. // Voltage Threshold 
#define Vpeak 20. // Spikes Peak 

#define m0 0.01
const double Tm[4] = {20.,10.,20.,20.} ; // Membrane time constants 
#define argIext 0 

#define IF_TRANSIENT_IEXT 0 
#define T_IEXT 1.0E3 
#define nbIext 1.0E3 

#define IF_JabLoop 0
#define Ax 0 
#define Ay 0 
#define IF_M0 0
#define Je0 3.0

//////////////////////////
// External Perturbation
//////////////////////////
#define IF_Prtr 1 
#define IF_WEAK 0 
#define PrtrPop 1
#define PrtrProfile 1

#define IF_PROP 0 
#define IF_WEAKPROP 0 

#define IF_PULSE 0 
#define PULSE 40 

#define IF_OPSIN 0 
#define SIGMA_EXT .1 
#define OpsPb 1. 

#define IF_TIMECOURSE 0 
#define Tc 1.E3 
#define dTc 1.2E3 
#define dPrtr .25 

//////////////////////////
// Interactions Structure
//////////////////////////
#define IF_RING 0 // cosine 
#define IF_GAUSS 0 // Gaussian 
#define IF_EXP 1 // Exponential 

#define L 3.0
#define DIM 1 // dimension of the ring 

#define IF_CstCrecCff 1
const double Sigma[4] ={.125,.075,.125,.075}; 
#define SigmaExt 0.10

#define IF_Dij 1 
const double Dij[16] ={1.0,1.5,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.,1.,1.,1.,1.,1.,1.};

#define IF_LIGHT_SCATTER 0
#define SIG_SCATTER L/8.0 

#define IF_SPEC 0 // with specific connections 
#define PHI0 0. 

//////////////////////////
#define IF_ADAPTATION 0 
#define Tad 100. 
#define Gad 1. 

//////////////////////////
#define IF_SHARED 0
#define PROP_SHARED .25
#define CLUSTER_SIZE .5

#endif
