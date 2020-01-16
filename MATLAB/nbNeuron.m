function [nbN N]= nbNeuron(nbpop,n,IF_Nk,p)

    nbN = zeros(nbpop,1) ;

    N = n*10000 ;

    for i=1:nbpop 
        nbN(i)=N./nbpop ;
    end

    if IF_Nk
        N = 76800 ;
        if(nbpop==2)            
            nbN(1)= 57600 ; 
            nbN(2)= 19200 ;             
        end
        if(nbpop==4)
            nbN(1)= 57600 ; 
            nbN(2)= 6400 ; 
            nbN(3)= 6400 ; 
            nbN(4)= 6400 ; 
        end 
    end

end                                     