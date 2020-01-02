#ifndef _NEURON_COUNTS 
#define _NEURON_COUNTS 

#define N_THREADS 512

#define nbpop 2 // number of populations 
#define N_NEURONS 120000ULL // 30000ULL//76800ULL //43200ULL // total number of neurons 
#define nbPref 10000.0 

#define popSize .5 // .3333333333333333333333333333333333333333 // .75 // proportion of neurons in excitatory pop 
#define IF_Nk 0 // 1 // different number of neurons in each pop then fix popSize 

#define IF_QUENCH 0
#define IF_PRES 0

#define IF_LARGE 0 
const char* AtoB = "IE" ; 

#define K 2000. 

#define IF_CHUNKS 0
#define NCHUNKS 4 
#define MAXNEURONS 19200ULL 
/* #define MAXNEURONS 15360ULL  */

#define IF_SHARED 0 
#define PROP_SHARED 0.0 
#define CLUSTER_SIZE 0.0 

#define IF_AUTA 0 
#define AUTA_Pop 0  
#define AUTA_Pb 0.0 
const double AUTA_p[4] = {1., 1., 0., 0.} ; 

// 4pop 
const double Sigma[4] = {0.125,.075,.125,.075} ; 
#define IF_Dij 1 
const double Dij[16] ={1.0,1.0,0.5,1.0,1.0,1.0,1.0,1.0,1.0,1.,1.,1.,1.,1.,1.,1.};

#define IF_LONG_RANGE 0 
#define AmpLR 1.0 
const double SigmaLR[4] = {.5,.5,.2250,.075} ; 
const double DijLR[16] = {1.0,0.0, 1.0,0.0, 1.0,1.0, 1.0,1.0,1.0 ,1.,1.,1., 1.,1.,1.,1.} ; 

#define L 3.0 // M_PI // size of the ring 
#define IF_RING 0 // standard ring with cosine interactions 
#define IF_SPACE 1 // Spatially structured connections 
#define IF_GAUSS 0 // Gaussian profile 
#define IF_EXP 1 // Exponential profile 

#define DIMENSION 1 // Dimension of the ring 
#define IF_PERIODIC 1 // periodic gaussian 

#define IF_SPEC 0 // sqrt(K) specific connections 

#define IF_MATRIX 0 // save Cij matrix 
#define IF_SPARSEVEC 1 // save sparse vectors 

#endif
