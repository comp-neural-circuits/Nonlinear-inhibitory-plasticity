%%%%%% The following code generates Fig. 6E,F from Miehl & Gjorgjieva 2022
%%%%%% PLoS CB. https://doi.org/10.1371/journal.pcbi.1010682

clear all

%% Parameter definitions

wEE=1; % initial E-to-E weight strength
wEI=1; % initial I-to-E weight strength
wIE=0.5; % initial E-to-I weight strength
cE=1; % E postsynaptic LTD/LTP threshold 
cI=1; % I postsynaptic LTD/LTP threshold 

NE=1; % Number of presynaptic E neurons
NI=1;

rhoE=2; % Presynaptic E rate in [Hz]
rhoI=0.5;  % External E rate onto I neurons in [Hz]
FR_E=max(NE*rhoE*wEE-NI*rhoI*wEI,0); % E postsynaptic firing rate in [Hz]
FR_I=rhoI+wIE*rhoE; % I firing rate in [Hz]

tau_FR_E=10; % Time constant for E neuron rate dynamics in [ms]
tau_FR_I=10; % Time constant for I neuron rate dynamics in [ms]
tau_wEE=1000; % Timescale for E plasticity in [ms]
tau_wEI=200; % Timescale for I plasticity in [ms]

total_time=20000; % total simulation time in [ms]

dt=0.1; % Integration timestep

counter=0;

save_timestep=1;

start_plot=5000;

%% Simulation start
for tt=dt:dt:total_time
    
    
    rhoE=0.5*sin(0.01*tt)+2;

    FR_E=FR_E+(-FR_E+max(NE*(rhoE)*wEE-NI*FR_I*wEI,0))/tau_FR_E*dt;
    
    FR_I=FR_I+(-FR_I+rhoI+wIE*rhoE)/tau_FR_I*dt;

    wEE=wEE+(rhoE*FR_E*(FR_E-cE))/tau_wEE*dt;
    wEI=wEI+(FR_I*FR_E*(FR_E-cI))/tau_wEI*dt;

    if mod(round(tt/dt)/10,save_timestep)==0
        counter=counter+1;
        save_stuff(counter,1)=FR_E;
        save_stuff(counter,2)=FR_I;
        save_stuff(counter,3)=wEE;
        save_stuff(counter,4)=wEI;
        save_stuff(counter,5)=rhoE;
    end 
end


%% Plot figures
map = brewermap(3,'Blues');
map2 = brewermap(3,'Reds');

width_of_lines=1;
size_font=8;

h2=figure;
subplot(6,3,1)
hold on
plot([dt:save_timestep:total_time-start_plot]./1000,save_stuff(start_plot+1:end,3),'Color',map(2,:),'LineWidth',width_of_lines)
plot([dt:save_timestep:total_time-start_plot]./1000,save_stuff(start_plot+1:end,4),'Color',map2(2,:),'LineWidth',width_of_lines)
hold off
XLABEL=xlabel('Time');
YLABEL=ylabel('w');


subplot(6,3,4)
hold on
plot([dt:save_timestep:total_time-start_plot]./1000,save_stuff(start_plot+1:end,1),'k','LineWidth',width_of_lines)
hold off
ylim([0,2])
XLABEL=xlabel('Time');
YLABEL=ylabel('Postsyn. rate');
set(gca,'linewidth',width_of_lines)



subplot(6,3,[2,5])

size_grid=5;
w1 = linspace(0,size_grid,20);
w2 = linspace(0,size_grid,20);

hold on
plot(save_stuff(start_plot+1:end,3),save_stuff(start_plot+1:end,4))
plot([0,size_grid],(NE.*max(save_stuff(:,5)).*[0,size_grid]-cE)./NI./(rhoI+wIE*max(save_stuff(:,5))),'k','Linewidth',width_of_lines)
plot([0,size_grid],(NE.*min(save_stuff(:,5)).*[0,size_grid]-cE)./NI./(rhoI+wIE*min(save_stuff(:,5))),'k','Linewidth',width_of_lines)
hold off
xlim([0,3])
ylim([0,3])
XLABEL=xlabel('w^{EE}');
YLABEL=ylabel('w^{EI}');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)

subplot(6,3,[3,6])

size_grid=5;
w1 = linspace(0,size_grid,20);
w2 = linspace(0,size_grid,20);

hold on
plot(save_stuff(start_plot+1:end,3),save_stuff(start_plot+1:end,4))
plot([0,size_grid],(NE.*max(save_stuff(:,5)).*[0,size_grid]-cE)./NI./(rhoI+wIE*max(save_stuff(:,5))),'k','Linewidth',width_of_lines)
plot([0,size_grid],(NE.*min(save_stuff(:,5)).*[0,size_grid]-cE)./NI./(rhoI+wIE*min(save_stuff(:,5))),'k','Linewidth',width_of_lines)
hold off
xlim([0.9,1.4])
ylim([0.5,1])
XLABEL=xlabel('w^{EE}');
YLABEL=ylabel('w^{EI}');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)