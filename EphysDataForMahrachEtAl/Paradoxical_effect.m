%
% by Guang Chen 2020.01.09
% PC: pyramidal neuron; PV: PV positive inhibitory interneuron
% ALM: anterior lateral motor cortex; S1: Somatosensory cortex
% data format in mat file
%
% --------------- key variables ---------------
% PC_rate_base: [Nunit x power x trial]
% single trial baseline spike rate of each unit under specific power of photo stimulation
%
% PC_rate_opto: [Nunit x power x trial]
% single trial spike rate of each unit during photo stimulation
% 
% PV_rate_base: [Nunit x power x trial]
% PV_rate_opto: [Nunit x power x trial]
%
% Power: [1 x N power level tested]  
% N power levels tested in the experiments, in mW 
% to compute light intensity mW/mm^2, the power is divided by the area of
% the light beam (e.g. for 1 mm diameter beam, light intensity is 'Power/(0.5^2*pi)')
% 10 different powers were tested in S1, but each neuron usually had 5 level of light stimulation
% 18 different powere were tested in ALM data, [1:9] for 1 mm beam diameter photostimulation [10:18] for 2 mm beam diameter photostimulation
% -----------------------------------------------------------------
%
% There are four mat files with data from different experiments:
% PC_PV_rate_power_ALM_superficial_20191231.mat
% PC_PV_rate_power_ALM_deep_20191231.mat
% PC_PV_rate_power_S1_20191231.mat
% PC_PV_rate_power0.5mwmm-2_S1_20191231.mat

%Run the following two code to reproduce the two experimental %figures in the paper
clear all
close all

%% Plot Figure 1 in the paper
clear;
path='.\';
figure;
for j=1:2% load ALM superficial/deep layer data
    if j==1
        % load ALM superficial layer data and plot figure1B
        PC_rate=[];
        PV_rate=[];
        filename='PC_PV_rate_power_ALM_superficial_20191231';
        load([path,filename,'.mat']);                       % PC_rate_base=[Nunit x power x trial] power dimension; rows 1:9 is 1mm diameter light spot sti.; rows 10:18 are 2mm sti.
        PC_rate_tmp=nanmean(PC_rate_base,3);                % trial average for base rate
        PC_rate(:,:,1)=repmat(nanmean(PC_rate_tmp,2),1,9);  % power condition average for base rate
        PC_rate(:,:,2)=nanmean(PC_rate_opto(:,10:end,:),3); % trial average for 2mm photo stimulation
        PV_rate_tmp=nanmean(PV_rate_base,3);                % trial average for base rate
        PV_rate(:,:,1)=repmat(nanmean(PV_rate_tmp,2),1,9);  % power condition average for base rate
        PV_rate(:,:,2)=nanmean(PV_rate_opto(:,10:end,:),3); % trial average for 2mm photo stimulation
    else
        % load ALM deep layer data and plot figure1C
        PC_rate=[];
        PV_rate=[];
        filename='PC_PV_rate_power_ALM_deep_20191231';
        load([path,filename,'.mat']);                       % PC_rate_base=[Nunit x power x trial] power dimension 1:9 are 1mm diameter light spot sti. 10:18 are 2mm sti.
        PC_rate_tmp=nanmean(PC_rate_base,3);                % trial average for base rate
        PC_rate(:,:,1)=repmat(nanmean(PC_rate_tmp,2),1,9);  % power condition average for base rate
        PC_rate(:,:,2)=nanmean(PC_rate_opto(:,10:end,:),3); % trial average for 2mm photo stimulation
        PV_rate_tmp=nanmean(PV_rate_base,3);                % trial average for base rate
        PV_rate(:,:,1)=repmat(nanmean(PV_rate_tmp,2),1,9);  % power condition average for base rate
        PV_rate(:,:,2)=nanmean(PV_rate_opto(:,10:end,:),3); % trial average for 2mm photo stimulation
    end
    
    % plot normalized rate (photo/baseline) for individual cells in figure 1B,C top
    subplot(3,3,j);
    hold on;
    plot(Power/(pi*1^2),(PC_rate(:,:,2)./PC_rate(:,:,1))','r');     % light intensity=power/area mW/mm2
    hold on;
    plot(Power/(pi*1^2),(PV_rate(:,:,2)./PV_rate(:,:,1))','b');
    hold on;
    plot([0.01 100],[1 1],'k--');hold on;
    set(gca, 'xscale', 'log');
    set(gca, 'yscale', 'log');
    xlim([0.01 100]);
    ylim([0.01 100]);
    text(10,10,'PC','color','r');hold on;
    text(10,50,'PV','color','b');hold on;
    if j==1
        title('ALM Layer 2/3');
    else
        title('ALM Layer 5');
    end
    
    % plot normalized mean rate across neurons and resample 10000 times in figure 1B,C bottom
    mean_norm_ratePC=[];
    for k=1:10000
        krand=randsample(size(PC_rate,1),size(PC_rate,1),1);
        mean_norm_ratePC=[mean_norm_ratePC;mean(squeeze(PC_rate(krand,:,2)))./mean(squeeze(PC_rate(krand,:,1)))];
    end
    mean_norm_ratePV=[];
    for k=1:10000
        krand=randsample(size(PV_rate,1),size(PV_rate,1),1);
        mean_norm_ratePV=[mean_norm_ratePV;mean(squeeze(PV_rate(krand,:,2)))./mean(squeeze(PV_rate(krand,:,1)))];
    end
    subplot(3,3,3+j);
    hold on;
    errorbar(Power/(pi*1^2),nanmean(mean_norm_ratePC)',nanstd(mean_norm_ratePC)','r','Linewidth',2);
    hold on;
    errorbar(Power/(pi*1^2),nanmean(mean_norm_ratePV)',nanstd(mean_norm_ratePV)','b','Linewidth',2);
    hold on;
    if j==1
        xlabel('Intensity(mW/mm^2)');
        ylabel('Relative rate photo/baseline');
    end
    plot([0.01 100],[1 1],'k--');hold on;
    set(gca, 'xscale', 'log');
    set(gca, 'yscale', 'log');
    xlim([0.01 100]);
    ylim([0.01 100]);
    
    % plot slope (relative spike rate change vs light intensity change) in
    % figure 1F (ALM Layer 5 pyramidal neuron and PV neuron)
    if j==2
        slopePC_tmp=(mean_norm_ratePC(:,3)-mean_norm_ratePC(:,1))/((Power(3)-Power(1))./(pi*1.^2));
        slopePV_tmp=(mean_norm_ratePV(:,3)-mean_norm_ratePV(:,1))/((Power(3)-Power(1))./(pi*1.^2));
        subplot(3,3,8);
        bar(1,mean(slopePC_tmp),'r');hold on;
        bar(2,mean(slopePV_tmp),'b');hold on;
        errorbar(1,mean(slopePC_tmp),std(slopePC_tmp),'r');hold on;
        errorbar(2,mean(slopePV_tmp),std(slopePV_tmp),'b');hold on;
        p=length(find((slopePC_tmp-slopePV_tmp)>0))/length(slopePC_tmp);
        if p>0.5
            p=1-p;
        end
        text(1,-1.5,['p=',num2str(p)]);
        xlabel('L5-PC     L5-PV');
        ylabel('slope');
        ylim([-2 1]);
        sloperatio=slopePV_tmp./slopePC_tmp;
        title(['ALM L5: slope ratio PV/PC=',num2str(mean(sloperatio)),'+-',num2str(std(sloperatio))]);
    end
end

% load ALM superficial and deep layer data and plot figure1E
obj1=load([path,'PC_PV_rate_power_ALM_superficial_20191231.mat']);
obj2=load([path,'PC_PV_rate_power_ALM_deep_20191231.mat']);
PV_rate_tmp=nanmean(obj1.PV_rate_base,3);                   % trial average for base rate
PV_rate1(:,:,1)=repmat(nanmean(PV_rate_tmp,2),1,9);         % power condition average for base rate
PV_rate1(:,:,2)=nanmean(obj1.PV_rate_opto(:,10:end,:),3);   % trial average for 2mm photo stimulation
PV_rate_tmp=nanmean(obj2.PV_rate_base,3);                   % trial average for base rate
PV_rate2(:,:,1)=repmat(nanmean(PV_rate_tmp,2),1,9);         % power condition average for base rate
PV_rate2(:,:,2)=nanmean(obj2.PV_rate_opto(:,10:end,:),3);   % trial average for 2mm photo stimulation

data1=PV_rate1(:,4,2)./PV_rate1(:,4,1);                     % relative spike rate (photo/baseline) at light intensity of 0.5mw/mm2
data2=PV_rate2(:,4,2)./PV_rate2(:,4,1);

subplot(3,3,7);
bar([1 2],[mean(data1) mean(data2)]);hold on;
errorbar(1,mean(data1),std(data1)/sqrt(length(data1)),'b');hold on; % mean+-sem
errorbar(2,mean(data2),std(data2)/sqrt(length(data2)),'b');hold on;
plot([0 3],[1 1],'k--');hold on;
[h,p]=ttest2(data1,data2);
text(1,2.5,['p=',num2str(p)]);
ylim([0 3]);
xlabel('L2/3     L5');
ylabel('Relative rate photo/baseline');
title('ALM');

% load S1 data and plot figure 1D,G
PC_rate=[];
PV_rate=[];
filename='PC_PV_rate_power_S1_20191231';
load([path,filename,'.mat']);                       % PC_rate_base=[Nunit*power*trial]
PC_rate(:,:,1)=nanmean(PC_rate_base,3);             % trial average for base rate
PC_rate(:,:,2)=nanmean(PC_rate_opto,3);             % trial average for photo stimulation rate
idx=find(sum(~isnan(squeeze(PC_rate(:,:,1))),2)>1); % find neurons with more than one level of photo stimulation
PC_rate=PC_rate(idx,:,:);
PV_rate(:,:,1)=nanmean(PV_rate_base,3);
PV_rate(:,:,2)=nanmean(PV_rate_opto,3);
idx=find(sum(~isnan(squeeze(PV_rate(:,:,1))),2)>1);
PV_rate=PV_rate(idx,:,:);

% plot normalized rate (photo/baseline) for individual cells in figure 1D top
subplot(3,3,3);
hold on;
for k=1:size(PC_rate,1)
    i_power=find(~isnan(PC_rate(k,:,2)));
    plot(Power(i_power)/(pi*1^2),(PC_rate(k,i_power,2)./PC_rate(k,i_power,1)),'r');%power/area mW/mm2
    hold on;
end
for k=1:size(PV_rate,1)
    i_power=find(~isnan(PV_rate(k,:,2)));
    plot(Power(i_power)/(pi*1^2),(PV_rate(k,i_power,2)./PV_rate(k,i_power,1)),'b');
    hold on;
end
plot([0.01 100],[1 1],'k--');hold on;
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
text(10,10,'PC','color','r');hold on;
text(10,50,'PV','color','b');hold on;
xlim([0.01 100]);
ylim([0.01 100]);
title('Somatosensory cortex');

% plot normalized mean rate across neurons and resample 10000 times in
% figure 1D bottom
mean_norm_ratePC=[];
for j=1:length(Power)
    idxtmp=find(~isnan(PC_rate(:,j,2)));
    PC_sti_tmp_j=PC_rate(idxtmp,j,2);
    PC_base_tmp_j=PC_rate(idxtmp,j,1);
    for k=1:10000
        krand=randsample(length(idxtmp),length(idxtmp),1);
        mean_norm_ratePC(k,j)=nanmean(PC_sti_tmp_j(krand))/nanmean(PC_base_tmp_j(krand));
    end
end
mean_norm_ratePV=[];
for j=1:length(Power)
    idxtmp=find(~isnan(PV_rate(:,j,2)));
    PV_sti_tmp_j=PV_rate(idxtmp,j,2);
    PV_base_tmp_j=PV_rate(idxtmp,j,1);
    for k=1:10000
        krand=randsample(length(idxtmp),length(idxtmp),1);
        mean_norm_ratePV(k,j)=nanmean(PV_sti_tmp_j(krand))/nanmean(PV_base_tmp_j(krand));
    end
end
subplot(3,3,6);
hold on;
errorbar(Power/(pi*1^2),nanmean(mean_norm_ratePC)',nanstd(mean_norm_ratePC)','r','Linewidth',2);
hold on;
errorbar(Power/(pi*1^2),nanmean(mean_norm_ratePV)',nanstd(mean_norm_ratePV)','b','Linewidth',2);
hold on;
plot([0.01 100],[1 1],'k--');hold on;
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
xlim([0.01 100]);
ylim([0.01 100]);

% plot slope in figure 1G
slopePC_tmp=[];
rate1=squeeze(nanmean(PC_rate(:,1:2,2),2));%opto one of the 1st and 2nd values is spike rate, the other is NaN
rate2=squeeze(nanmean(PC_rate(:,1:2,1),2));%base
rate3=squeeze(nanmean(PC_rate(:,3:4,2),2));%opto
rate4=squeeze(nanmean(PC_rate(:,3:4,1),2));%base
for k=1:10000
    krand=randsample(length(rate1),length(rate1),1);
    slopePC_tmp=[slopePC_tmp (nanmean(rate3(krand))/nanmean(rate4(krand))-nanmean(rate1(krand))/nanmean(rate2(krand)))/((mean(Power(3:4))-mean(Power(1:2)))/(pi*1^2))];
end
slopePV_tmp=[];
rate1=squeeze(nanmean(PV_rate(:,1:2,2),2));%opto
rate2=squeeze(nanmean(PV_rate(:,1:2,1),2));%base
rate3=squeeze(nanmean(PV_rate(:,3:4,2),2));%opto
rate4=squeeze(nanmean(PV_rate(:,3:4,1),2));%base
for k=1:10000
    krand=randsample(length(rate1),length(rate1),1);
    slopePV_tmp=[slopePV_tmp (nanmean(rate3(krand))/nanmean(rate4(krand))-nanmean(rate1(krand))/nanmean(rate2(krand)))/((mean(Power(3:4))-mean(Power(1:2)))/(pi*1^2))];
end
subplot(3,3,9);
bar(1,mean(slopePC_tmp),'r');hold on;
bar(2,mean(slopePV_tmp),'b');hold on;
errorbar(1,mean(slopePC_tmp),std(slopePC_tmp),'r');hold on;
errorbar(2,mean(slopePV_tmp),std(slopePV_tmp),'b');hold on;
p=length(find((slopePC_tmp-slopePV_tmp)>0))/length(slopePC_tmp);
if p>0.5
    p=1-p;
end
text(1,-1.8,['p=',num2str(p)]);
xlabel('S1-PC     S1-PV');
ylabel('slope');
ylim([-2 1]);
sloperatio=slopePV_tmp./slopePC_tmp;
title(['S1: slope ratio PV/PC=',num2str(mean(sloperatio)),'+-',num2str(std(sloperatio))]);



%% Plot Figure 2 in the paper
clear all;
path='.\';
figure;
for j=1:3       %for ALM superficial/deep/S1 data
    for z=1:2   %for PC/PV
        switch j
            case 1
                % load and plot ALM superficial layer data
                PC_rate=[];
                PV_rate=[];
                filename='PC_PV_rate_power_ALM_superficial_20191231';
                load([path,filename,'.mat']);                       % PC_rate_base=[Nunit x power x trial] power dimension 1:9 are 1mm sti., 10:18 is 2mm sti
                PC_rate_tmp=nanmean(PC_rate_base,3);                % trial average for base rate
                PC_rate(:,:,1)=repmat(nanmean(PC_rate_tmp,2),1,9);  % power condition average for base rate
                PC_rate(:,:,2)=nanmean(PC_rate_opto(:,10:end,:),3); % trial average for 2mm photo stimulation
                PV_rate_tmp=nanmean(PV_rate_base,3);                % trial average for base rate
                PV_rate(:,:,1)=repmat(nanmean(PV_rate_tmp,2),1,9);  % power condition average for base rate
                PV_rate(:,:,2)=nanmean(PV_rate_opto(:,10:end,:),3); % trial average for 2mm photo stimulation
                if z==1
                    xdata=squeeze(PC_rate(:,5,1));                  % take baseline spike rate at light intensity of 0.5 mw/mm2
                    ydata=squeeze(PC_rate(:,5,2));                  % take photostimulation spike rate at light intensity of 0.5 mw/mm2
                else
                    xdata=squeeze(PV_rate(:,5,1));                  % take baseline spike rate at light intensity of 0.5 mw/mm2
                    ydata=squeeze(PV_rate(:,5,2));                  % take photostimulation spike rate at light intensity of 0.5 mw/mm2
                end
            case 2
                % load and plot ALM deep layer data
                PC_rate=[];
                PV_rate=[];
                filename='PC_PV_rate_power_ALM_deep_20191231';
                load([path,filename,'.mat']);                       % PC_rate_base=[Nunit x power x trial] power dimension [1:9] is 1mm sti. [10:18] is 2mm sti
                PC_rate_tmp=nanmean(PC_rate_base,3);                % trial average for base rate
                PC_rate(:,:,1)=repmat(nanmean(PC_rate_tmp,2),1,9);  % power condition average for base rate
                PC_rate(:,:,2)=nanmean(PC_rate_opto(:,10:end,:),3); % trial average for 2mm opto stimulation
                PV_rate_tmp=nanmean(PV_rate_base,3);                % trial average for base rate
                PV_rate(:,:,1)=repmat(nanmean(PV_rate_tmp,2),1,9);  % power condition average for base rate
                PV_rate(:,:,2)=nanmean(PV_rate_opto(:,10:end,:),3); % trial average for 2mm opto stimulation
                if z==1
                    xdata=squeeze(PC_rate(:,5,1));                  % take baseline spike rate at light intensity of 0.5 mw/mm2
                    ydata=squeeze(PC_rate(:,5,2));                  % take photostimulation spike rate at light intensity of 0.5 mw/mm2
                else
                    xdata=squeeze(PV_rate(:,5,1));                  % take baseline spike rate at light intensity of 0.5 mw/mm2
                    ydata=squeeze(PV_rate(:,5,2));                  % take photostimulation spike rate at light intensity of 0.5 mw/mm2
                end
            case 3
                % load and plot S1 data
                filename='PC_PV_rate_power0.5mwmm-2_S1_20191231';
                load([path,filename,'.mat']);           %PC_rate_base=[Nunit x power x trial] light intensity=0.5mw/mm2
                if z==1
                    xdata=nanmean(PC_rate_base,2);      %take baseline spike rate at light intensity of 0.5 mw/mm2
                    idx=find(~isnan(xdata));
                    xdata=xdata(idx);
                    ydata=nanmean(PC_rate_opto,2);      %take photostimulation spike rate at light intensity of 0.5 mw/mm2
                    ydata=ydata(idx);
                else
                    xdata=nanmean(PV_rate_base,2);
                    idx=find(~isnan(xdata));
                    xdata=xdata(idx);
                    ydata=nanmean(PV_rate_opto,2);
                    ydata=ydata(idx);
                end
        end
        
        xdata(find(xdata<=0.01))=0.01;
        xdata(find(xdata>=100))=100;
        ydata(find(ydata<=0.01))=0.01;
        ydata(find(ydata>=100))=100;
        
        % scatter plot, spike rate photostimulation vs spike rate baseline
        subplot(4,3,j+(z-1)*6);
        if z==1
            scatter(xdata,ydata,'r');hold on;%for PC
        else
            scatter(xdata,ydata,'b');hold on;%for PV
        end
        plot([0.01 100],[0.01 100],'k--');hold on;
        set(gca, 'xscale', 'log');
        set(gca, 'yscale', 'log');
        xlabel('Spike rate, baseline(spikes/s)','FontSize',7);
        ylabel('Spike rate, photostimulation','FontSize',7);
        if j==1&z==1
            title('ALM Layer 2/3, PC');
        end
        if j==2&z==1
            title('ALM Layer 5, PC');
        end
        if j==3&z==1
            title('S1, PC');
        end
        if z==2
            title('PV');
        end
        
        % pie plot
        for k=1:1000
            krand=randsample(length(xdata),length(xdata),1);
            group1=find(ydata(krand)<0.1);                                                      % group 1: rate <0.1 under pv photoinhibition
            group2=intersect(find(ydata(krand)>=0.1),find(abs(ydata(krand)-xdata(krand))<0.1)); % group 2: spike rate no change, abs rate change <0.1 Hz under photoinhibition
            group3=intersect(find(ydata(krand)>=0.1),find((ydata(krand)-xdata(krand))>=0.1));   % group 3: spike rate increase >=0.1 Hz under photoinhibition
            group4=intersect(find(ydata(krand)>=0.1),find((ydata(krand)-xdata(krand))<=-0.1));  % group 4: spike rate decrease <=-0.1 Hz under photoinhibition
            frac1(k)=length(group1)/length(krand);
            frac2(k)=length(group2)/length(krand);
            frac3(k)=length(group3)/length(krand);
            frac4(k)=length(group4)/length(krand);
        end
        subplot(4,3,j+(z-1)*6+3);
        pie([mean(frac1) mean(frac2) mean(frac3) 1-(mean(frac1)+mean(frac2)+mean(frac3))]);
        
        % display proportion for each of the groups. 
        clear *_label
        if z==1 
            neuron_population_label = 'PC population, ';
        elseif z==2
            neuron_population_label = 'PV population, ';
        end
        if j==1
            area_label = 'ALM L2/3, ';
        elseif j==2
            area_label = 'ALM L5, ';
        elseif j==3
            area_label = 'S1, ';
        end
        disp(['=========== Figure 2 pie chart ',neuron_population_label, area_label, ' ==============='])
        disp(['Neuron% with spike rate <0.1 Hz upon PV stim: ',num2str(round(mean(frac1)*1000)/10),'%+-',num2str(round(std(frac1)*1000)/10),'%']);
        disp(['Neuron% with no change in spike rate upon PV stim: ',num2str(round(mean(frac2)*1000)/10),'%+-',num2str(round(std(frac2)*1000)/10),'%']);
        disp(['Neuron% with spike rate increase (>=0.1Hz) upon PV stim: ',num2str(round(mean(frac3)*1000)/10),'%+-',num2str(round(std(frac3)*1000)/10),'%']);
        disp(['Neuron% with spike rate decrease (<=0.1Hz) upon PV stim: ',num2str(round((1-(mean(frac1)+mean(frac2)+mean(frac3)))*1000)/10),'%+-',num2str(round(std(frac4)*1000)/10),'%']);
            
        
    end
end

subplot(4,3,6);
legend('group 1','group 2','group 3','group 4')






