%%%%%% The following code generates Fig. 3D-F from Miehl and Gjorgjieva 2021


close all
clear all


%% Parameter definitions

wEE=1.5; % initial E-to-E weight strength
wEI=0.5; % initial I-to-E weight strength
wIE=0.5; % initial E-to-I weight strength
wEE2=1.5; % initial E-to-E weight strength
wEI2=0.5; % initial I-to-E weight strength
cE0=1; 
cI0=1;
cE=cE0*0.7; % initial E postsynaptic LTD/LTP threshold
cI=cI0*1.3; % initial i postsynaptic LTD/LTP threshold
cE2=cE0;
cI2=cI0;

NE=1; % Number of presynaptic E neurons
NI=1;

rhoE=2; % Presynaptic E rate in [Hz]
rhoE2=2;
rhoI=0.5; % External E rate onto I neurons in [Hz]
FR_E=max(NE*rhoE*wEE-NI*rhoI*wEI,0); % E postsynaptic firing rate in [Hz]
FR_I=rhoI+wIE*rhoE; % I firing rate in [Hz]
FR_E2=max(NE*rhoE*wEE2-NI*rhoI*wEI2,0); 
FR_I2=rhoI+wIE*rhoE; 

tau_FR_E=10; % Time constant for E neuron rate dynamics in [ms]
tau_FR_I=10; % Time constant for I neuron rate dynamics in [ms]
tau_wEE=1000; % Timescale for E plasticity in [ms]
tau_wEI=200; % Timescale for I plasticity in [ms]
tau_c=2; % Timescale threshold dynamics in [ms]

total_time=15000; % total simulation time in [ms]

dt=0.1; % Integration timestep

counter=0;
vec_rhoE=[0:0.01:5];

%% Simulation start
for tt=dt:dt:total_time
    
    
    counter=counter+1;

    FR_E=FR_E+(-FR_E+max(NE*rhoE*wEE-NI*FR_I*wEI,0))/tau_FR_E*dt;
    FR_I=FR_I+(-FR_I+rhoI+wIE*rhoE)/tau_FR_I*dt;

    Delta_wEE=(rhoE*FR_E*(FR_E-cE))/tau_wEE*dt;
    Delta_wEI=(FR_I*FR_E*(FR_E-cI))/tau_wEI*dt;
    wEE=wEE+Delta_wEE;
    wEI=wEI+Delta_wEI;
    wEE(wEE<0)=0;
    wEI(wEI<0)=0;

    cE=cE+Delta_wEE/tau_c*dt;
    cI=cI-Delta_wEI/tau_c*dt;
    
    save_stuff(counter,1)=FR_E;
    save_stuff(counter,2)=FR_I;
    save_stuff(counter,3)=wEE;
    save_stuff(counter,4)=wEI;
    save_stuff(counter,5)=cE;
    save_stuff(counter,6)=cI;
    
    FR_E2=FR_E2+(-FR_E2+max(NE*rhoE2*wEE2-NI*FR_I2*wEI2,0))/tau_FR_E*dt;
    FR_I2=FR_I2+(-FR_I2+rhoI+wIE*rhoE2)/tau_FR_I*dt;

    Delta_wEE2=(rhoE2*FR_E2*(FR_E2-cE2))/tau_wEE*dt;
    Delta_wEI2=(FR_I2*FR_E2*(FR_E2-cI2))/tau_wEI*dt;
    wEE2=wEE2+Delta_wEE2;
    wEI2=wEI2+Delta_wEI2;
    wEE2(wEE2<0)=0;
    wEI2(wEI2<0)=0;
    
    cE2=cE2+Delta_wEE2/tau_c*dt;
    cI2=cI2-Delta_wEI2/tau_c*dt;
    

    save_stuff2(counter,1)=FR_E2;
    save_stuff2(counter,2)=FR_I2;
    save_stuff2(counter,3)=wEE2;
    save_stuff2(counter,4)=wEI2;
    save_stuff2(counter,5)=cE2;
    save_stuff2(counter,6)=cI2;
end


map = brewermap(3,'Blues');
map0 = brewermap(4,'Set1');
map2 = brewermap(3,'Reds');


width_of_lines=1;
size_font=8;

%% Plot figures
h1=figure;

subplot(3,3,1)
x_ax_postFR=[0:0.01:2];
hold on
    plot(x_ax_postFR,(rhoE.*x_ax_postFR.*(x_ax_postFR-cE0*0.7))./tau_wEE,'Color',map(2,:),'LineWidth',width_of_lines);
    plot(x_ax_postFR,((rhoI+wIE*rhoE).*x_ax_postFR.*(x_ax_postFR-cI0*1.3))./tau_wEI,'Color',map2(2,:),'LineWidth',width_of_lines);
    plot(x_ax_postFR,zeros(length(x_ax_postFR),1),'black','LineWidth',width_of_lines)
    plot([cE0*0.7,cE0*0.7],[-1,2.5],':k','LineWidth',width_of_lines)
    plot([cI0*1.3,cI0*1.3],[-1,2.5],':k','LineWidth',width_of_lines)
hold off
ylim([-4,4.5]*10^(-3))
XLABEL=xlabel('Postsyn. rate \nu^E in [Hz]');
YLABEL=ylabel('\Delta w');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
hLegend=legend({'\Delta w^{EE}','\Delta w^{EI}',},'FontSize',8,'FontName','Arial','location','northwest');
hLegend.ItemTokenSize = [15,18];
legend('boxoff')
set(gca,'linewidth',width_of_lines)

subplot(3,3,2)
hold on
plot([dt:dt:total_time]./1000,save_stuff(:,5),'--','Color',map(2,:),'LineWidth',width_of_lines)
plot([dt:dt:total_time]./1000,save_stuff(:,6),'--','Color',map2(2,:),'LineWidth',width_of_lines)
plot([dt:dt:total_time]./1000,save_stuff2(:,5),'Color',map(2,:),'LineWidth',width_of_lines)
plot([dt:dt:total_time]./1000,save_stuff2(:,6),'Color',map2(2,:),'LineWidth',width_of_lines)
hold off
ylim([0.7,1.3])
XLABEL=xlabel('Time in [s]');
YLABEL=ylabel('Threshold c');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)

subplot(3,3,3)
yyaxis left
hold on
plot([dt:dt:total_time]./1000,save_stuff(:,3),'--','Color',map(2,:),'LineWidth',width_of_lines)
plot([dt:dt:total_time]./1000,save_stuff(:,4),'--','Color',map2(2,:),'LineWidth',width_of_lines)
plot([dt:dt:total_time]./1000,save_stuff2(:,3),'Color',map(2,:),'LineWidth',width_of_lines)
plot([dt:dt:total_time]./1000,save_stuff2(:,4),'Color',map2(2,:),'LineWidth',width_of_lines)
hold off
ylim([0,8])
XLABEL=xlabel('Time in [s]');
YLABEL=ylabel('w');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)

yyaxis right
hold on
plot([dt:dt:total_time]./1000,save_stuff(:,1),'--k','LineWidth',width_of_lines)
plot([dt:dt:total_time]./1000,save_stuff2(:,1),'k','LineWidth',width_of_lines)
ylim([0.75,1.25])
YLABEL=ylabel('Postsynaptic rate in [Hz]');




h1.Renderer='Painters';
%print(h1,'dynamic_thresholds_c_V3','-dpdf','-bestfit')