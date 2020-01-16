% Routine used to plot scatter plots and pie charts from the
% simulation data

clear all ; 
GlobalVars 

baseline = ImportData(model, nbpop, dir,'IdvRates', N, K, g, IF_IEXT, prtrPop, Iext ) ; 

try
    for i=1:length(baseline(1,:))-1
        IdvRates(i) = mean(baseline(:,i+1)) ;
    end
catch
    fprintf('Error Rates not found\n')
    return ;
end

prtr = ImportData(model, nbpop, dir,'IdvRates', N, K, g, IF_IEXT, prtrPop, Iprtr);

try
    for i=1:length(prtr(1,:))-1
        IdvRatesPrtr(i) = mean(prtr(:,i+1)) ;
    end
catch
    fprintf('Error prtr Rates not found\n')
    return ;
end

for i=1:2
                
    MeanRate(i) = mean( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ;
    VarRate(i) = var( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ; 

    MeanRatePrtr(i) = mean( IdvRatesPrtr( Cpt(i)+1:Cpt(i+1) ) ) ;
    VarRatePrtr(i) = var( IdvRatesPrtr( Cpt(i)+1:Cpt(i+1) ) ) ; 
    
    Baseline = IdvRates( Cpt(i)+1: Cpt(i+1) ); 
    RatesPrtr = IdvRatesPrtr( Cpt(i)+1: Cpt(i+1) ) ; 

    ROI1 = find(Baseline<=.01) ; 
    Baseline(ROI1)=.01 ; 
    ROI2 = find(RatesPrtr<=.01) ; 
    RatesPrtr(ROI2)=.01 ; 

    ROI = intersect(ROI1,ROI2) ; 
    ZeroIdx = union(ROI,ROI2) ; 

    SupIdx = ( RatesPrtr-Baseline ) ; %./ ( Baseline + RatesPrtr ) ; 
    Ratio = SupIdx ./ Baseline ;
    IdxZo = [] ;

    if(~isempty(ZeroIdx))
        IdxZo = ZeroIdx ;
    else
        IdxZo = [] ;
    end

    if(any(abs(Ratio)<=0.01))
        IdxNc = setdiff(find(abs(Ratio)<=0.01),IdxZo) ; 
    else
        IdxNc = [] ; 
    end
    
    if(any(SupIdx>0)) 
        IdxUp = setdiff(setdiff(find(SupIdx>0),IdxZo),IdxNc) ; 
    else
        IdxUp = [] ; 
    end
    if(any(SupIdx<0))
        IdxDn = setdiff(setdiff(find(SupIdx<0),IdxZo),IdxNc) ; 
    else
        IdxDn = [] ;
    end

    propSum=0 ;
    if(~isempty(SupIdx)) 
        propUp = round( length(IdxUp) / length(SupIdx) * 100 ) ;
        propDn = round( length(IdxDn) / length(SupIdx) * 100 ) ; 
        propZo = round( length(IdxZo) / length(SupIdx) * 100 ) ;
        propNc = round( length(IdxNc) / length(SupIdx) * 100 ) ; 
        propSum = propUp + propDn + propNc + propZo ;        
    else
        propUp = 0 ;
        propZo = 100 ;
        propDn = 0 ;
        propNc = 0 ;
        propSum = propUp + propDn + propNc + propZo ;
    end

    fprintf('%s propUp %.0f propDn %.0f propZo %.0f propNc %.0f sum %.0f\n', popList(i) , propUp, propDn, propZo, propNc, propSum ) 

    while(propSum~=100)
        if(propSum>100)
            if(propZo>=1)
                propZo = propZo - 1 ;
            elseif(propZo<1)
                if(propNc>=1)
                    propNc = propNc -1 ;
                end
            end
        else 
            if(propZo>=1)
                propZo = propZo + 1 ;
            elseif(propZo<1)
                if(propNc>=1)
                    propNc = propNc + 1 ;
                end
            end
        end

        propSum = propUp + propDn + propNc + propZo ;
    end

    fprintf('%s propUp %.0f propDn %.0f propZo %.0f propNc %.0f sum %.0f\n', popList(i) , propUp, propDn, propZo, propNc, propSum ) 


    Idx = randi([Cpt(i)+1 Cpt(i+1)],300,1) ; 
    Baseline = IdvRates(Idx) ; 
    RatesPrtr = IdvRatesPrtr(Idx) ; 

    ROI1 = find(Baseline<=.01) ; 
    Baseline(ROI1)=.01 ; 
    ROI2 = find(RatesPrtr<=.01) ; 
    RatesPrtr(ROI2)=.01 ; 

    ROI = intersect(ROI1,ROI2) ; 
    ZeroIdx = union(ROI,ROI2) ; 

    SupIdx = ( RatesPrtr-Baseline ) ; %./ ( Baseline + RatesPrtr ) ; 
    Ratio = SupIdx ./ Baseline ;
    
    IdxZo = [] ;

    if(~isempty(ZeroIdx))
        IdxZo = ZeroIdx ;
    else
        IdxZo = [] ;
    end

    if(any(abs(Ratio)<=0.01))
        IdxNc = setdiff(find(abs(Ratio)<=0.01),IdxZo) ; 
    else
        IdxNc = [] ; 
    end
    
    if(any(SupIdx>0)) 
        IdxUp = setdiff(setdiff(find(SupIdx>0),IdxZo),IdxNc) ; 
    else
        IdxUp = [] ; 
    end
    if(any(SupIdx<0))
        IdxDn = setdiff(setdiff(find(SupIdx<0),IdxZo),IdxNc) ; 
    else
        IdxDn = [] ;
    end
 
    figtitle=sprintf('Scatter%s_Iprtr%.2f',popList(i),Iprtr(prtrPop)) ; 

    fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ;
    xlabel('Baseline (Hz)')
    ylabel('Light On (Hz)')
    % end

    scatter(Baseline,RatesPrtr,20,cl{i},'filled','MarkerFaceAlpha',.5,'LineWidth',.5) ; 
    
    plot([.01 100],[.01 100],'--k','Linewidth',.5) 
    xlim([.01 100])
    ylim([.01 100])
    set(gca,'yscale','log')
    set(gca,'xscale','log')

    drawnow ;
    hold off ; 

    figtitle=sprintf('PieProp%s_Iprtr%.3f',popList(i),Iprtr(prtrPop)) ; 
    fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ;
    X = [propUp, propNc, propDn, propZo] ;
    labels = {'Up','NC','Down','Zero'} ;
    labels = {'','','',''} ; 
    %pie(X,labels) 
    p = pie(X) ; 
    pieColor = gray(4) ; 
    isProp = logical([propUp propNc propDn propZo]) ;
    for l=4:-1:1
        if(~isProp(l))
            pieColor(l,:) = [] ;
        end
    end
    colormap(pieColor) 

    axis off ; 
    drawnow ; 
    hold off ; 
end

fprintf('Norm. Rates: ')
fprintf('%.3f | ', MeanRatePrtr./MeanRate)
fprintf('\n')
