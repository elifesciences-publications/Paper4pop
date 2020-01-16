% Global variables to update before using other routine.
% needs to be supplemented with the simulation data from ../../LIF/Simulations/*
warning off ; 

%%%%%%%%%%%%%%%%%%%%%

model = 'LIF' ; % network model
nbpop = 2 ; % number of populations 2 or 4.
dir = '' % name of dir with set of parameters : '2pop', 'model1', 'model11', 'model2'.
file = 'MeanRates' ; % see CheckBal

%%%%%%%%%%%%%%%%%%%%% 

N = 7 ; % first digit in number of neurons 
K = 500 ; % average number of inputs per neuron and per population
g=1 ; % gain prefactor 
IF_Nk = 1 ; % if 1 different number of neurons in each population

J = ImportJab(model,nbpop,dir) ; % Connectivity parameters 
Iext = ExternalInput(model,nbpop,dir) ; % external input
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ; % number of neurons in each population
Cpt = CptNeuron(nbpop,nbN) ; % cumulative number of neurons in each population

%%%%%%%%%%%%%%%%%%%%%

mkSize = 5 ; % markersize
mk = 'o' ; % markertype
alp = 1 ; % markerFaceAlpha
popList = ['E' 'I' 'S' 'V'] ;
cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ; 

%%%%%%%%%%%%%%%%%%%%%

THRESHOLD = .1 ; % rate threshold 
CV_THRESHOLD = .1 ; % CV threshold

%%%%%%%%%%%%%%%%%%%%%

IF_NORM = 1 ; % Normalizes rates to baseline
FIGPERPOP = 1 ; % Separates PC/PV from SOM/VIP
IF_IDVTraces = 1 ; % adds idindividual traces
nbIdv = [10 10 10 10] ; % number of individual traces 

IF_POWER = 1 ; % Power profile 
I0 = 8. ; 
P0 = .5 ; 

IDX = 2 ; % Skips first data points
IF_LOGSCALE = 1 ; % log-log plot
IF_LOGSCALEX = 0 ; % log scale on x axis only

%%%%%%%%%%%%%%%%%%%%%

IF_IEXT = 'Delta' ; % Profile of the perturbation 
prtrPop = 2 ; % Idx of the perturbed population 
Iprtr = Iext ; % Total external input (feedforward + perturbation)

prtrAmp = 0 ; % Amplitude of the perturbation; set to .28 and .45 for the scatter plots and piecharts 

if(prtrPop>0)
    Iprtr(prtrPop) = Iprtr(prtrPop) + prtrAmp ; 
else
    Iprtr = prtrAmp * ones(1,nbpop) ;
end

v_Iprtr = 0:.1:1. ; % Range of the perturbation intensity

%%%%%%%%%%%%%%%%%%%%%
