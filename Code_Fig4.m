%%%%%% The following code generates Fig. 4 from Miehl and Gjorgjieva 2021



close all
clear all

%% Parameter definitions

wEE=1.5; % initial E-to-E weight strength
wEI=0.5; % initial I-to-E weight strength
wIE=0.5; % initial E-to-I weight strength
wEE2=1.5; % initial E-to-E weight strength
wEI2=0.5; % initial I-to-E weight strength
cE=1; % E postsynaptic LTD/LTP threshold 
cI=1; % I postsynaptic LTD/LTP threshold 

NE=1; % Number of presynaptic E neurons
NI=1;

rhoE=2; % Presynaptic E rate in [Hz]
rhoE2=2;
rhoE_0=2;
rhoI=0.5; % External E rate onto I neurons in [Hz]
FR_E=max(NE*rhoE*wEE-NI*rhoI*wEI,0); % E postsynaptic firing rate in [Hz]
FR_I=rhoI+wIE*rhoE; % I firing rate in [Hz]
FR_E2=max(NE*rhoE*wEE2-NI*rhoI*wEI2,0); 
FR_I2=rhoI+wIE*rhoE; 


tau_FR_E=10; % Time constant for E neuron rate dynamics in [ms]
tau_FR_I=10; % Time constant for I neuron rate dynamics in [ms]
tau_wEE=1000; % Timescale for E plasticity in [ms]
tau_wEI=200; % Timescale for I plasticity in [ms]

total_time=2000; % total simulation time in [ms]
time_change=1000; % timepoint of perturbation in [ms]

dt=0.1; % Integration timestep

counter=0;

vec_rhoE=[0:0.01:5];


%% Simulation start
for tt=dt:dt:total_time
    
        if round(tt,1)==(time_change+dt)

        for jj=1:length(vec_rhoE)
            rhoE_loop=vec_rhoE(jj);
            nuI_loop=vec_rhoE(jj);
            
            % delta wE dependent on rhoE
            delta_wE_case5_2=rhoE_loop*max(NE*rhoE_loop*wEE-(rhoI+rhoE_loop*wIE)*wEI,0)*(max(NE*rhoE_loop*wEE-(rhoI+rhoE_loop*wIE)*wEI,0)-cE)/tau_wEE;
            save_wE_case5(jj,1)=delta_wE_case5_2;
            
            
            % delta wI dependent on nuI
            %delta_wI_case5_2=rhoI_loop*max(NE*rhoE*wEE-(rhoI_loop+rhoE*wIE)*wEI,0)*(max(NE*rhoE*wEE-(rhoI_loop+rhoE*wIE)*wEI,0)-cI)/tau_wEI;
            delta_wI_case5_2=nuI_loop*max(NE*rhoE*wEE-(nuI_loop)*wEI,0)*(max(NE*rhoE*wEE-(nuI_loop)*wEI,0)-cI)/tau_wEI;
            save_wI_case5(jj,1)=delta_wI_case5_2;
            
         
        end      
        rhoE=1.5; % here the E input rate is perturbed
        rhoE2=2.5;
    end
    
    counter=counter+1;

    FR_E=FR_E+(-FR_E+max(NE*rhoE*wEE-NI*FR_I*wEI,0))/tau_FR_E*dt;
    FR_I=FR_I+(-FR_I+rhoI+wIE*rhoE)/tau_FR_I*dt;

    wEE=wEE+(rhoE*FR_E*(FR_E-cE))/tau_wEE*dt;
    wEI=wEI+(FR_I*FR_E*(FR_E-cI))/tau_wEI*dt; % linear I plasticity rule
    wEE(wEE<0)=0;
    wEI(wEI<0)=0;
    
    
    FR_E2=FR_E2+(-FR_E2+max(NE*rhoE2*wEE2-NI*FR_I2*wEI2,0))/tau_FR_E*dt;
    FR_I2=FR_I2+(-FR_I2+rhoI+wIE*rhoE2)/tau_FR_I*dt;
    
    wEE2=wEE2+(rhoE2*FR_E2*(FR_E2-cE))/tau_wEE*dt;
    wEI2=wEI2+(FR_I2*FR_E2*(FR_E2-cI))/tau_wEI*dt; % linear I plasticity rule
    wEE2(wEE2<0)=0;
    wEI2(wEI2<0)=0;

    save_stuff(counter,1)=FR_E;
    save_stuff(counter,2)=FR_I;
    save_stuff(counter,3)=wEE;
    save_stuff(counter,4)=wEI;
    save_stuff2(counter,1)=FR_E2;
    save_stuff2(counter,2)=FR_I2;
    save_stuff2(counter,3)=wEE2;
    save_stuff2(counter,4)=wEI2;
    
end

for jj=1:length(vec_rhoE)
    rhoE_loop=vec_rhoE(jj);
    nuI_loop=vec_rhoE(jj);

    % delta wE dependent on rhoE
    delta_wE_case5_2=rhoE_loop*max(NE*rhoE_loop*wEE-(rhoI+rhoE_loop*wIE)*wEI,0)*(max(NE*rhoE_loop*wEE-(rhoI+rhoE_loop*wIE)*wEI,0)-cE)/tau_wEE;
    save_wE_case5(jj,2)=delta_wE_case5_2;

    % delta wI dependent on nuI
    %delta_wI_case5_2=rhoI_loop*max(NE*rhoE*wEE-(rhoI_loop+rhoE*wIE)*wEI,0)*(max(NE*rhoE*wEE-(rhoI_loop+rhoE*wIE)*wEI,0)-cI)/tau_wEI;
    delta_wI_case5_2=nuI_loop*max(NE*rhoE*wEE-(nuI_loop)*wEI,0)*(max(NE*rhoE*wEE-(nuI_loop)*wEI,0)-cI)/tau_wEI;
    save_wI_case5(jj,2)=delta_wI_case5_2;
    
    % delta wE dependent on rhoE 2
    delta_wE_case5_22=rhoE_loop*max(NE*rhoE_loop*wEE2-(rhoI+rhoE_loop*wIE)*wEI2,0)*(max(NE*rhoE_loop*wEE2-(rhoI+rhoE_loop*wIE)*wEI2,0)-cE)/tau_wEE;
    save_wE_case52(jj,2)=delta_wE_case5_22;



    % delta wI dependent on rhoI 2
    %delta_wI_case5_22=rhoI_loop*max(NE*rhoE2*wEE2-(rhoI_loop+rhoE2*wIE2)*wEI2,0)*(max(NE*rhoE2*wEE2-(rhoI_loop+rhoE2*wIE2)*wEI2,0)-cI)/tau_wEI;
    delta_wI_case5_22=nuI_loop*max(NE*rhoE2*wEE2-(nuI_loop)*wEI2,0)*(max(NE*rhoE2*wEE2-(nuI_loop)*wEI2,0)-cI)/tau_wEI;
    save_wI_case52(jj,2)=delta_wI_case5_22;

end      


%% Plot figures
map = brewermap(3,'Blues');
map0 = brewermap(4,'Set1');
map2 = brewermap(3,'Reds');
map3 = brewermap(6,'Greens');

width_of_lines=1;
size_font=8;

upperlim=2/tau_wEE;
lowerlim=-1.1/tau_wEE;

h2=figure;


subplot(3,3,2)
hold on
plot([500+dt:dt:total_time]./1000,save_stuff(5001:end,1),'--','Color',map(2,:),'LineWidth',width_of_lines)
plot([500+dt:dt:total_time]./1000,save_stuff(5001:end,2),'--','Color',map2(2,:),'LineWidth',width_of_lines)
plot([500+dt:dt:total_time]./1000,save_stuff2(5001:end,1),'Color',map(2,:),'LineWidth',width_of_lines)
plot([500+dt:dt:total_time]./1000,save_stuff2(5001:end,2),'Color',map2(2,:),'LineWidth',width_of_lines)
plot([time_change,time_change]./1000,[0,3.1],':k','LineWidth',width_of_lines)
hold off
ylim([0,2.5])
XLABEL=xlabel('Time in [s]');
YLABEL=ylabel('Firing rate in [Hz]');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
hLegend=legend({'E decr.','I decr.','E incr.','I incr.'},'FontSize',size_font,'FontName','Arial','location','southeast');
hLegend.ItemTokenSize = [15,18];
legend('boxoff')
set(gca,'linewidth',width_of_lines)

subplot(3,3,5)
hold on
plot([500+dt:dt:total_time]./1000,save_stuff(5001:end,3),'--','Color',map(2,:),'LineWidth',width_of_lines)
plot([500+dt:dt:total_time]./1000,save_stuff(5001:end,4),'--','Color',map2(2,:),'LineWidth',width_of_lines)
plot([500+dt:dt:total_time]./1000,save_stuff2(5001:end,3),'Color',map(2,:),'LineWidth',width_of_lines)
plot([500+dt:dt:total_time]./1000,save_stuff2(5001:end,4),'Color',map2(2,:),'LineWidth',width_of_lines)
plot([time_change,time_change]./1000,[0,2.5],':k','LineWidth',width_of_lines)
hold off
ylim([0,2.5])
%axis square
XLABEL=xlabel('Time in [s]');
YLABEL=ylabel('w');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
%set(h3,'Units','Inches');
%set(gcf, 'PaperPositionMode','auto');
%hLegend=legend('w^{EE} decr.','w^{EI} decr.','w^{EE} incr.','w^{EI} incr.','location', 'North' );
%legend('boxoff')
set(gca,'linewidth',width_of_lines)
%pos2 = get(h3,'Position');
%set(gca, 'box', 'off')
%set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos2(3), pos2(4)]);
  

subplot(3,3,3)
hold on
plot(vec_rhoE,save_wE_case5(:,2),'Color',map(1,:),'LineWidth',width_of_lines)
plot(vec_rhoE,save_wE_case5(:,1),'Color',map(2,:),'LineWidth',width_of_lines)
plot(vec_rhoE,save_wE_case52(:,2),'Color',map(3,:),'LineWidth',width_of_lines)
plot(vec_rhoE,zeros(jj,1),'black','LineWidth',width_of_lines)
hold off
xlim([0 3.5])
ylim([lowerlim,upperlim])
%axis square
XLABEL=xlabel('Presyn. exc. rate \rho^E in [Hz]');
YLABEL=ylabel('\Delta w^{EE}');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
%set(h1,'Units','Inches');
%set(gcf, 'PaperPositionMode', 'auto');
hLegend=legend({'\rho_{disr}^E=1.5','\rho_{base}^E=2','\rho_{disr}^E=2.5'},'FontSize',size_font,'FontName','Arial','location','northwest');
hLegend.ItemTokenSize = [15,18];
legend('boxoff')
set(gca,'linewidth',width_of_lines)
% pos2 = get(h1,'Position');
% set(gca, 'box', 'off')
% set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos2(3), pos2(4)]);  

subplot(3,3,6)
hold on
plot(vec_rhoE,save_wI_case5(:,2),'Color',map2(1,:),'LineWidth',width_of_lines)
plot(vec_rhoE,save_wI_case5(:,1),'Color',map2(2,:),'LineWidth',width_of_lines)
plot(vec_rhoE,save_wI_case52(:,2),'Color',map2(3,:),'LineWidth',width_of_lines)
plot(vec_rhoE,zeros(jj,1),'black','LineWidth',width_of_lines)
hold off
xlim([0 3.5])
ylim([-0.005,0.035])
%axis square
XLABEL=xlabel('Presyn. inh. rate \nu^I in [Hz]');
YLABEL=ylabel('\Delta w^{EI}');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
%set(h4,'Units','Inches');
%set(gcf, 'PaperPositionMode', 'auto');
hLegend=legend({'\rho_{disr}^E=1.5','\rho_{base}^E=2','\rho_{disr}^E=2.5'},'FontSize',size_font,'FontName','Arial','location','northeast');
hLegend.ItemTokenSize = [15,18];
legend('boxoff')
set(gca,'linewidth',width_of_lines)


subplot(3,3,4)
hold on
plot([0,5],(NE*rhoE.*[0,5]-cE)./NI./(rhoI+wIE*rhoE),'color',[0.8,0.8,0.8],'Linewidth',width_of_lines)
plot([0,5],(NE*rhoE_0.*[0,5]-cE)./NI./(rhoI+wIE*rhoE_0),'color',[0.4,0.4,0.4],'Linewidth',width_of_lines)
plot([0,5],(NE*rhoE2.*[0,5]-cE)./NI./(rhoI+wIE*rhoE2),'k','Linewidth',width_of_lines)
hold off
XLABEL=xlabel('w^{EE}');
YLABEL=ylabel('w^{EI}');
xlim([0,5])
ylim([0,5])
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
hLegend=legend({'\rho_{disr}^E=1.5','\rho_{base}^E=2','\rho_{disr}^E=2.5'},'FontSize',size_font,'FontName','Arial','location','northwest');
hLegend.ItemTokenSize = [15,18];
legend('boxoff')
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)



