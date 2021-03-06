%clear all ; 
GlobalVars 

Iext = ExternalInput(model,nbpop,dir) ;    
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

Cth=1000 ;

try
    data = ImportData(model,nbpop,dir,'Raster',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iprtr) ; 
    Spikes = sortrows(data) ; 
    Spikes(:,2) = Spikes(:,2)./1000. ;
catch
    fprintf('FILE NOT FOUND\n') ;
    return ;
end

WinSize = ceil( Spikes(end,2) ) ; 
fprintf('Duration %ds\n', WinSize ) 

reSizeIdx = 1:length( unique( Spikes ) ) ; 

[nbSpk nidx]= hist(Spikes(:,1),unique( Spikes(:,1) ) ) ; 
CumSumSpk = [0 cumsum(nbSpk)] ; 

CVpop = NaN(Cpt(nbpop+1),1) ;

fprintf('nbN ')
for j=1:nbpop
    fprintf('%d ', nbN(j))
end 

fprintf('\nlength ROI ') 
for j=1:nbpop 
    X=[] ; 
    Y=[] ; 
    for i=1:nbN(j) 
        X(i) = L * mod( double(i), sqrt( double( nbN(j) ) ) ) / sqrt( double( nbN(j) ) ) ; 
        Y(i) = L * floor( double(i) / sqrt( double( nbN(j) ) ) ) / sqrt( double( nbN(j) ) ) ; 
    end 
    ROI{j} = find( (X-L/2).^2 + (Y-L/2).^2 <= Cth.^2 / 4 ) ; 
    fprintf('%d ',length(ROI{j}))
end
fprintf('\n')

SpkTimes = {} ; 
for i=1:length(CumSumSpk)-1 
    fprintf('# %d nIdx %d nbSpk %.3f ', nidx(i), reSizeIdx(i), nbSpk(i)) 
    try 
        SpkTimes = [ SpkTimes;{Spikes(CumSumSpk(i)+1:CumSumSpk(i+1),2).'} ] ; 
        % SpkTimes = [ SpkTimes;{Spikes(CumSumSpk(i)+1:CumSumSpk(i)+nbSpk(i),2).'} ] ; 
        ISI = SpkTimes{i}(2:nbSpk(i))-SpkTimes{i}(1:nbSpk(i)-1) ; 
        % ISI = SpkTimes{i}(2:nbSpk(i))-SpkTimes{i}(1:nbSpk(i)-1) ; 
        CV(i) = sqrt(var(ISI)) ./ mean(ISI) ; 
    catch 
        CV(i) = nan ; 
    end
    fprintf('CV %.3f ',CV(i))
    fprintf('\r')

    for j=1:nbpop
        if( CV(i)>CV_THRESHOLD && ismember(nidx(i)-Cpt(j),ROI{j}) ) 
            CVpop( nidx(i) ) = CV(i) ; 
        end
    end
end
fprintf('\n')

fprintf('mean CV ') 

for i=1:nbpop
    switch i 
      case 1
        figname='CVdistE' ;
      case 2
        figname='CVdistI' ;
      case 3
        figname='CVdistS' ;
      case 4
        figname='CVdistV' ;
    end

    CVx = CVpop( Cpt(i)+1:Cpt(i+1) ) ; 
    CVx(isnan(CVx)) = [] ;
    fprintf('%.2f ', mean( CVx ) ) 
    
    % if( ishandle( findobj('type','figure','name',figname) ) ) 
    %     fig = findobj('type','figure','name',figname) ; 
    %     figure(fig); hold on ; 
    % else
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        xlabel('CV') 
        ylabel('pdf') 
    % end
    
    h = histogram( CVx , 30, 'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i},'EdgeAlpha',alp,'Linewidth',2) ; 
    xlim([0 2]) 
    drawnow ; 
    
    if(IF_SAVE) 
        figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_IEXT) ;
        figdir= sprintf('%s/Baseline',figdir) ;
        fprintf('Writing %s \n',figdir) 
        try
            mkdir(figdir)
        end
        
        ProcessFigure(fig, fullfile(figdir,figname), 2.2, [1.33*2.2, 2.2]) ;
    end        

    % hold off ;
end
fprintf('\n')
    
% try
%     data = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iext(prtrPop)) ;
% catch
%     return ;
% end

% for i=1:length(data(1,:))-1
%     IdvRates(i) = mean(data(:,i+1)) ;
% end

% for i=1:nbpop
    
%     X = IdvRates(Cpt(i)+1:Cpt(i+1) ) ; 
%     Y = CVpop(Cpt(i)+1:Cpt(i+1) ) ; 
    
%     IDX = ~isnan(Y) ; 
%     Y = Y(IDX) ; 
%     X = X(IDX) ;     
%     fprintf('%.2f ', mean(Y) ) ;
    
%     figname = sprintf('RatesVsCV%s',popList(i)) ; 
%     fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
%     xlabel('Rates') 
%     ylabel('CV') 
%     plot(X,Y,'+','color',cl{i}) 
% end
% fprintf('\n')