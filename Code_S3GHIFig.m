%%%%%% The following code generates S3G-I Fig. from Miehl & Gjorgjieva 2022
%%%%%% PLoS CB. https://doi.org/10.1371/journal.pcbi.1010682

clear all

wEE=0.7; % initial E-to-E weight strength
wEI=0.5; % initial I-to-E weight strength
wIE=0.5; % initial E-to-I weight strength
wEE2=0.7; % initial E-to-E weight strength
wEI2=0.5; % initial I-to-E weight strength
wIE2=0.5; % initial E-to-I weight strength
cE=1; % E postsynaptic LTD/LTP threshold 
cI=1; % I postsynaptic LTD/LTP threshold 

NE=1; % Number of presynaptic E neurons
NI=1;
NE2=1; % Number of presynaptic E neurons
NI2=1;

rhoE=2; % Presynaptic E rate in [Hz]
rhoE2=2;
rhoI=1; % External E rate onto I neurons in [Hz]
rhoI2=1; % External E rate onto I neurons in [Hz]
FR_I=rhoI+wIE*rhoE; % I firing rate in [Hz]
FR_I2=rhoI2+wIE2*rhoE2; 
FR_E=max(NE*rhoE*wEE+NE2*rhoE2*wEE2-NI*FR_I*wEI-NI2*FR_I2*wEI2,0); % E postsynaptic firing rate in [Hz]

%FR_E2=max(NE*rhoE*wEE2-NI*rhoI*wEI2,0); 


tau_FR_E=10; % Time constant for E neuron rate dynamics in [ms]
tau_FR_I=10; % Time constant for I neuron rate dynamics in [ms]
tau_wEE=1000; % Timescale for E plasticity in [ms]
tau_wEI=200; % Timescale for I plasticity in [ms]
tau_c=2; % Timescale threshold dynamics in [ms]

total_time=40000; % total simulation time in [ms]

dt=0.1; % Integration timestep

counter=0;
vec_rhoE=[0:0.01:8];

%% Simulation start
for tt=dt:dt:total_time
    
    if round(tt,1)==(20000+dt)
        wEIbalancebeginning=wEE/wEI;
        for jj=1:length(vec_rhoE)
            rhoE_loop=vec_rhoE(jj);
            rhoI_loop=vec_rhoE(jj);
            
            % delta wE dependent on rhoE
            delta_wE_case5_2=rhoE_loop*max(NE*rhoE_loop*wEE+NE2*rhoE2*wEE2-NI*(rhoI+rhoE_loop*wIE)*wEI-NI2*(rhoI2+rhoE2*wIE2)*wEI2,0)*(max(NE*rhoE_loop*wEE+NE2*rhoE2*wEE2-NI*(rhoI+rhoE_loop*wIE)*wEI-NI2*(rhoI2+rhoE2*wIE2)*wEI2,0)-cE)/tau_wEE;
            save_wE_case5(jj,1)=delta_wE_case5_2;
            
             % delta wE dependent on rhoI
             delta_wE_case5_3=rhoE*max(NE*rhoE*wEE+NE2*rhoE2*wEE2-NI*(rhoI_loop+rhoE*wIE)*wEI-NI2*(rhoI_loop+rhoE2*wIE2)*wEI2,0)*(max(NE*rhoE*wEE+NE2*rhoE2*wEE2-NI*(rhoI_loop+rhoE*wIE)*wEI-NI2*(rhoI_loop+rhoE2*wIE2)*wEI2,0)-cE)/tau_wEE;
             save_wE_case5_3(jj,1)=delta_wE_case5_3;
            
        end      
        rhoE=1; % Perturb E input
        rhoE2=3.5;
    end
    
    counter=counter+1;

    FR_E=FR_E+(-FR_E+max(NE*rhoE*wEE+NE2*rhoE2*wEE2-NI*FR_I*wEI-NI2*FR_I2*wEI2,0))/tau_FR_E*dt;
    FR_I=FR_I+(-FR_I+rhoI+wIE*rhoE)/tau_FR_I*dt;
    FR_I2=FR_I2+(-FR_I2+rhoI2+wIE2*rhoE2)/tau_FR_I*dt;

    Delta_wEE=(rhoE*FR_E*(FR_E-cE))/tau_wEE*dt;
    Delta_wEI=(FR_I*FR_E*(FR_E-cI))/tau_wEI*dt;
    
    wEE=wEE+Delta_wEE;
    wEI=wEI+Delta_wEI;
    wEE(wEE<0)=0;
    wEI(wEI<0)=0;
    
    Delta_wEE2=(rhoE2*FR_E*(FR_E-cE))/tau_wEE*dt;
    Delta_wEI2=(FR_I2*FR_E*(FR_E-cI))/tau_wEI*dt;
    
    wEE2=wEE2+Delta_wEE2;
    wEI2=wEI2+Delta_wEI2;
    wEE2(wEE2<0)=0;
    wEI2(wEI2<0)=0;
    
    cE=cE+(Delta_wEE+Delta_wEE2)/tau_c*dt;
    cI=cI-(Delta_wEI+Delta_wEI2)/tau_c*dt;
    

    save_stuff(counter,1)=FR_E;
    save_stuff(counter,2)=FR_I;
    save_stuff(counter,3)=wEE;
    save_stuff(counter,4)=wEI;
    save_stuff(counter,5)=cE;
    save_stuff(counter,6)=cI;
    
    save_stuff2(counter,1)=FR_E;
    save_stuff2(counter,2)=FR_I2;
    save_stuff2(counter,3)=wEE2;
    save_stuff2(counter,4)=wEI2;
    

end

for jj=1:length(vec_rhoE)
    rhoE_loop=vec_rhoE(jj);
    rhoI_loop=vec_rhoE(jj);

    delta_wE_case5_2=rhoE_loop*max(NE*rhoE_loop*wEE+NE2*rhoE2*wEE2-NI*(rhoI+rhoE_loop*wIE)*wEI-NI2*(rhoI2+rhoE2*wIE2)*wEI2,0)*(max(NE*rhoE_loop*wEE+NE2*rhoE2*wEE2-NI*(rhoI+rhoE_loop*wIE)*wEI-NI2*(rhoI2+rhoE2*wIE2)*wEI2,0)-cE)/tau_wEE;
    save_wE_case5(jj,2)=delta_wE_case5_2;

    delta_wE_case5_22=rhoE_loop*max(NE*rhoE*wEE+NE2*rhoE_loop*wEE2-NI*(rhoI+rhoE*wIE)*wEI-NI2*(rhoI2+rhoE_loop*wIE2)*wEI2,0)*(max(NE*rhoE*wEE+NE2*rhoE_loop*wEE2-NI*(rhoI+rhoE*wIE)*wEI-NI2*(rhoI2+rhoE_loop*wIE2)*wEI2,0)-cE)/tau_wEE;
    save_wE_case52(jj,2)=delta_wE_case5_22;
    
end      

wEIbalanceInput1after=wEE/wEI;
wEIbalanceInput2after=wEE2/wEI2;

%% Plot figures

map = brewermap(3,'Blues');
map2 = brewermap(3,'Reds');
h1=figure;

width_of_lines=1;
size_font=8;

subplot(3,3,1)
hold on
plot(vec_rhoE,save_wE_case5(:,1),'Color',map(2,:),'LineWidth',width_of_lines)
plot(vec_rhoE,save_wE_case5(:,2),'Color',map(3,:),'LineWidth',width_of_lines)
plot(vec_rhoE,save_wE_case52(:,2),'Color',map(1,:),'LineWidth',width_of_lines)
plot(vec_rhoE,zeros(length(vec_rhoE),1),'black','LineWidth',width_of_lines)
hold off
xlim([0 8])
ylim([-10^(-3),3*10^(-3)])
XLABEL=xlabel('Presyn. exc. rate \rho^E');
YLABEL=ylabel('\Delta w^{EE}');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',8);
set(gca,'FontSize',8,'FontName','Arial');
hLegend=legend({'Baseline','Input 1','Input 2'},'FontSize',8,'FontName','Arial','location','northwest');
hLegend.ItemTokenSize = [15,18];
legend('boxoff')
set(gca,'linewidth',width_of_lines)

subplot(3,3,2)
hold on
plot([18000:dt:total_time]./1000,save_stuff(180000:end,5),'--','Color',map(2,:),'LineWidth',width_of_lines)
plot([18000:dt:total_time]./1000,save_stuff(180000:end,6),'--','Color',map2(2,:),'LineWidth',width_of_lines)
plot([20000,20000]./1000,[0,2.5],':k','LineWidth',width_of_lines)
hold off
ylim([0.98,1.01])
xlim([18 total_time/1000])
XLABEL=xlabel('Time');
YLABEL=ylabel('c');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)

subplot(3,3,3)
hold on
plot([18000:dt:total_time]./1000,save_stuff(180000:end,3),'--','Color',map(2,:),'LineWidth',width_of_lines)
plot([18000:dt:total_time]./1000,save_stuff(180000:end,4),'--','Color',map2(2,:),'LineWidth',width_of_lines)
plot([18000:dt:total_time]./1000,save_stuff2(180000:end,3),'Color',map(2,:),'LineWidth',width_of_lines)
plot([18000:dt:total_time]./1000,save_stuff2(180000:end,4),'Color',map2(2,:),'LineWidth',width_of_lines)
plot([20000,20000]./1000,[0,2.5],':k','LineWidth',width_of_lines)
hold off
ylim([0.4,0.8])
xlim([18 total_time/1000])
XLABEL=xlabel('Time');
YLABEL=ylabel('w');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)
