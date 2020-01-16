function data = ImportData(model,nbpop,dir,file,n,k,g,IF_Prtr,nPrtr,Xprtr)

    if(nargin<12) 
        nPrtr = 2 ;
        if(nargin<10)
            IF_Prtr = '' ;
            Xprtr = [] ;
        end
    end

    popList = ['E' 'I' 'S' 'V'] ;
    if(nbpop==1)
        popList = 'I' ;
    end
    
    if strcmp(model,'Connectivity') 
        path = sprintf(['../%s/%dpop/N%d/K%d'],model,nbpop,n,k) ;
    else
        path = sprintf(['../%s/Simulations/%dpop/%s/N%d/K%d/g%.2f'],model,nbpop,dir,n,k,g) ; 
    end

    switch IF_Prtr
        
      case {'Delta'}
        if(nPrtr>0)
            add_path = sprintf(['%sPrtr/Iext_%s%.4f'], IF_Prtr, popList(nPrtr), Xprtr(nPrtr)) ;
        else
            add_path = sprintf(['%sPrtr/Iext_All%.4f'], IF_Prtr, Xprtr(1)) ;
        end

      case 'TIMECOURSE'
        Iext = ExternalInput(model,nbpop,dir) ; 
        Xprtr1 = Xprtr(nPrtr) ; 
        add_path = sprintf(['DeltaPrtr/Iext_%s%.4f/TIMECOURSE'], popList(nPrtr), Xprtr1) ; 

      case 'JabLoop' 
        add_path = sprintf(['DeltaPrtr/Iext_%s%.4f/Je0%.4fJee%.4f'], ...
                           popList(nPrtr), Xprtr(1), Xprtr(2), Xprtr(3)) ;
         
      otherwise
        add_path = '' ;
    end
    
    path = sprintf(['%s/%s'],path,add_path) ;     
    if strcmp(model,'Connectivity') 
        file = sprintf(['%s/%s.dat'],path,file) ;
        data = file ; 
    else
        file = sprintf(['%s/%s.txt'],path,file) ;

        fprintf('Reading %s \n',file) 
        
        data = [] ;
        try
            data = importdata(file,separator) ; 
        catch
            fprintf('DATA NOT FOUND\n') 
        end
    end
    
end