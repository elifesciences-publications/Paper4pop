warning off ; 

%%%%%%%%%%%%%%%%%%%%%

popList = ['E' 'I' 'S' 'V'] ;
cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ; 

%%%%%%%%%%%%%%%%%%%%%

model = 'LIF' ; 
nbpop = 2 ; 
dir = 'L5' ; 
file = 'MeanRates' ; 

%%%%%%%%%%%%%%%%%%%%% 

N = 3 ; 
g = 1 ; 
IF_Nk = 0 ; 

J = ImportJab(model,nbpop,dir) ; 
Iext = ExternalInput(model,nbpop,dir) ;    
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

nbIdv = [10 10 10 10] ; 

%%%%%%%%%%%%%%%%%%%%%

IF_LOOP = 0 ;
if(~IF_LOOP) 
    K = 500 ; 
    mkSize = 5 ;     
    IF_ROBUST = 0 ; 
    mk = 'o' ; 
    alp = 1 ; 
end

mk = 'o' ; 

%%%%%%%%%%%%%%%%%%%%%

THRESHOLD = .1 ; 
CV_THRESHOLD = .1 ; 

%%%%%%%%%%%%%%%%%%%%%

IF_NORM = 1 ; 
FIGPERPOP = 0 ; 
IF_IDVTraces = 1 ; 
IF_COUNTPIF = 0 ; 

IF_POWER = 2 ; 
I0 = 8. ; 
P0 = .5 ; 

IDX = 2 ; 
IF_CORRECTION = 0 ; 
IF_LOGSCALE = 1 ; 
IF_LOGSCALEX = 0 ; 

IF_MF_RATES = 0 ; 
IF_SAVE = 1 ; 

%%%%%%%%%%%%%%%%%%%%%

IF_PROP = 0 ; 
IF_PROPWEAK = 0 ; 

%%%%%%%%%%%%%%%%%%%%%

IF_IEXT = 'Gauss' ; 
prtrPop = 2 ; 
Iprtr = Iext ; 

prtrAmp = 0 ; %.5 ; % .28 et .45

if(prtrPop>0)
    Iprtr(prtrPop) = Iprtr(prtrPop) + prtrAmp ; 
else
    Iprtr = prtrAmp * ones(1,nbpop) ;
end

v_Iprtr = 0:.1:1. ; 
%v_Iprtr = [.1:.1:.9,1:1:10 ]; 

%%%%%%%%%%%%%%%%%%%%%

IF_RING = 'Exp' ; 
L = 3 ; 
DIM = 1 ; 

IF_SPACELOOP = exist('IF_SPACELOOP') ;
if(~IF_SPACELOOP)
    Cff = [.1] ; 
end

Crec = [.125 .075 0.125 .075] ; 

IF_Dij = 0 ; 
Dij = [1.0 1.0 1.6667 1.0] ; 

Cth = 100 ; 

% Cff = [] ; 
% Crec = [.25 .25 .25 .075] ; 
% Cth = 100 ; 

v_Cff = .075:.025:.375 ;  
%v_Cff = .2:.025:.375 ; 

%%%%%%%%%%%%%%%%%%%%%
