% Routine used to plot Norm. activities vs Gamma opto from the
% simulation data 
clear all ; 
GlobalVars

for i=1:length(v_Iprtr) 
    data = ImportData(model, nbpop, dir, 'IdvRates', N, K, g, IF_IEXT, prtrPop, Iext + v_Iprtr(i) ) ; 
    try
        for j=1:length(data(1,:))-1 
            IdvRates(i,j) = mean(data(:,j+1)) ;
        end
        
        fprintf(' Rates ') 
        for j=1:nbpop 
            m(j,i) = 0 ; 
            Rates = IdvRates(i, Cpt(j)+1:Cpt(j+1)) ; 
            if(i==1) 
                BL = IdvRates(i, Cpt(j)+1:Cpt(j+1)) ; 
            end
            m(j,i) = mean(Rates(find(BL>=THRESHOLD))) ; 
            fprintf('%.3f ', round(m(j,i),3) ) 
        end 
        fprintf('\n') 

    catch
        for j=1:nbN(nbpop)
            IdvRates(i,j) = nan ;
        end

        fprintf(' Rates ') 
        for j=1:nbpop 
            m(j,i) = nan ;
            Rates = IdvRates(i, Cpt(j)+1:Cpt(j+1)) ;
            fprintf('%.3f ', m(j,i) )
        end 
        fprintf('\n')        
    end
end

if(~FIGPERPOP)
    figtitle = sprintf('%s_RatesVsIopto',dir) ; 
    if(IF_LOGSCALE) 
        figtitle = sprintf('%s_RatesVsIopto_Log',dir) ;
    else
        figtitle = sprintf('%s_RatesVsIopto_LogX',dir) ;
    end
    if(IF_POWER) 
        figtitle = sprintf('%s_%d',figtitle,IF_POWER) ; 
    end
    if( ishandle( findobj('type','figure','name',figtitle) ) )
        fig = findobj('type','figure','name',figtitle) ; 
        fig = figure(fig); hold on ; 
    else
        fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
        xlabel('I_{opto}') 
        if(IF_NORM)
            ylabel('Norm. Rates') 
        else
            ylabel('Rates') 
        end
        if(IF_POWER) 
            xlabel('\Gamma_{opto} (mW/mm^2)') 
        end
    end
end

if IF_POWER
    v_Iprtr = P0 .* ( exp( v_Iprtr .* 20 ./ I0 ) - 1 ) ; 
end

for i=1:nbpop    
    if(FIGPERPOP)
        if(i==1 || i==2) 
            figtitle = sprintf('%s_RatesVsIopto%s_EI',dir,popList(prtrPop)) ; 
            if(IF_LOGSCALE)
                if(IF_IDVTraces)
                    figtitle = sprintf('%s_RatesVsIopto%sIdvLogEI',dir,popList(prtrPop)) ;
                else
                    figtitle = sprintf('%s_RatesVsIopto%sLogEI',dir,popList(prtrPop)) ;                    
                end
            end
        else 
            figtitle = sprintf('%s_RatesVsIopto%s_SV',dir,popList(prtrPop)) ; 
            if(IF_LOGSCALE)
                if(IF_IDVTraces)
                    figtitle = sprintf('%s_RatesVsIopto%sIdvLogSV',dir,popList(prtrPop)) ;
                else
                    figtitle = sprintf('%s_RatesVsIopto%sLogSV',dir,popList(prtrPop)) ;
                end
            end
        end
        if(IF_POWER) 
            figtitle = sprintf('%s_%d',figtitle,IF_POWER) ; 
            xlabel('\Gamma_{opto} (mW/mm^2)') 
        end
        fig = popPerFig(i,dir,figtitle) ;
        xlabel('I_{opto}') 
        if(IF_NORM)
            ylabel('Norm. Rates') 
        else 
            ylabel('Rates') 
        end 
        if(IF_POWER) 
            xlabel('\Gamma_{opto} (mW/mm^2)') 
        end
    end

    if(IF_NORM)
        NormRates = m(i,:)./m(i,1) ;
    else
        NormRates = m(i,:) ;
    end

    NormRates(find(NormRates<.01)) = .001 ;     
    id = find( ~isnan(NormRates) ) ; 
    id = IDX:length(v_Iprtr) ; 
    
    patchline(abs(v_Iprtr(id)), NormRates(id), 'linestyle','-','edgecolor',cl{i},'edgealpha',alp,'linewidth',1.5) 
    plot(abs(v_Iprtr(id)), NormRates(id), mk,'MarkerEdgeColor',cl{i},'MarkerSize',mkSize,'MarkerFaceColor','none','LineWidth', alp) 
    

    if( (i==2 || i==4 ) && IF_NORM) 
        if(IF_MF_RATES)
            if(i==2) 
                for j=1:2 
                    MFpop = MFrates(j,:)./MFrates(j,1) ;
                    plot(abs(v_Iprtr(id)), MFpop(id), '--','Color',cl{j}) 
                end
            elseif(i==4)
                for j=3:4
                    MFpop = MFrates(j,:)./MFrates(j,1) ;
                    plot(abs(v_Iprtr(id)), MFpop(id), '--','Color',cl{j}) 
                end
            end

        end
        
        plot(abs(v_Iprtr(IDX:end)),ones(1, length(v_Iprtr(IDX:end)) ),'--','Color','k') 
        
    end
    
    if(IF_IDVTraces) 

        countMax = 0 ;
        while countMax<nbIdv(i)
            nId = randi([Cpt(i)+1 Cpt(i+1)]) ; 
            if IdvRates(1,nId)>=THRESHOLD 
                countMax = countMax + 1 ; 
                IdvNormRates = IdvRates(:,nId)./IdvRates(1,nId) ;
                IdvNormRates(find(IdvNormRates<.01)) = .001 ;    
                patchline(v_Iprtr(id), IdvNormRates(id), 'linestyle','-','edgecolor',cl{i},'edgealpha',.2,'linewidth',1.5) 
            end 
        end 
        
    end
    
    if(IF_LOGSCALE) 
        xlim([.01 100])
        set(gca,'Xscale', 'log')
        if(prtrPop==1)
            ylim([.1 10])
        else
            ylim([.01 10])
        end
        set(gca,'Yscale', 'log') 
    else
        if(IF_LOGSCALEX)
            if(nbpop==2)
                xlim([.01 2])
            else
                xlim([.01 100])
            end
            set(gca,'Xscale', 'log')
        end
        if(i==2) 
            ylim([0 3]) 
        elseif(i==4) 
            ylim([0 4]) 
        end
    end    
end

hold off ;
