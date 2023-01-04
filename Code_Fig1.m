%%%%%% The following code generates Fig. 1 from Miehl & Gjorgjieva 2022
%%%%%% PLoS CB. https://doi.org/10.1371/journal.pcbi.1010682

close all
clear all

%% Parameter definitions

wEE=1.5; % initial E-to-E weight strength
wEI=0.5; % initial I-to-E weight strength
wIE=0.5; % initial E-to-I weight strength
wEE2=2.5; % initial E-to-E weight strength [unstable case]
wEI2=1; % initial I-to-E weight strength [unstable case]
cE=1; % E postsynaptic LTD/LTP threshold 
cI=1; % I postsynaptic LTD/LTP threshold 

NE=1; % Number of presynaptic E neurons
NI=1;

rhoE=2; % Presynaptic E rate in [Hz]
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

dt=0.1; % Integration timestep

counter=0;
vec_rhoE=[0:0.01:5]; % vector of E input firing rates in [Hz]

%% Calculate Fig. 1C
 vec_nuI=[0:1:6]; % vector of I firing rates in [Hz]
% for bb2=1:length(vec_nuI)
%     for bb=1:length(vec_rhoE) 
%         
%         FR_E_loop=max(NE*vec_rhoE(bb)*wEE-vec_nuI(bb2)*wEI,0);
%         save_dwEE(bb2,bb)=(vec_rhoE(bb)*FR_E_loop*(FR_E_loop-cE))/tau_wEE;
%         
%     end
% end

vec_wEI=[0:0.25:1.5]; % vector of wEI
for bb2=1:length(vec_wEI)
    for bb=1:length(vec_rhoE) 
        FR_I_loop=rhoI+wIE*vec_rhoE(bb);
        FR_E_loop=max(NE*vec_rhoE(bb)*wEE-NI*FR_I_loop*vec_wEI(bb2),0);
        save_dwEE(bb2,bb)=(vec_rhoE(bb)*FR_E_loop*(FR_E_loop-cE))/tau_wEE;
        
    end
end

%% Simulation start
for tt=dt:dt:total_time
    
    counter=counter+1;

    FR_E=FR_E+(-FR_E+max(NE*rhoE*wEE-NI*FR_I*wEI,0))/tau_FR_E*dt;
    FR_I=FR_I+(-FR_I+rhoI+wIE*rhoE)/tau_FR_I*dt;

    wEE=wEE+(rhoE*FR_E*(FR_E-cE))/tau_wEE*dt;
    wEI=wEI+(FR_I*(FR_E-cI))/tau_wEI*dt; % linear I plasticity rule
    wEE(wEE<0)=0;
    wEI(wEI<0)=0;
    
    
    FR_E2=FR_E2+(-FR_E2+max(NE*rhoE*wEE2-NI*FR_I2*wEI2,0))/tau_FR_E*dt;
    FR_I2=FR_I2+(-FR_I2+rhoI+wIE*rhoE)/tau_FR_I*dt;
    
    wEE2=wEE2+(rhoE*FR_E2*(FR_E2-cE))/tau_wEE*dt;
    wEI2=wEI2+(FR_I2*(FR_E2-cI))/tau_wEI*dt; % linear I plasticity rule
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


%% Plot figures
map = brewermap(3,'Blues');
map2 = brewermap(3,'Reds');
map3 = brewermap(6,'Greens');
map4 = brewermap(length(vec_nuI),'Blues');

width_of_lines=1;
size_font=8;

h1=figure;
subplot(3,3,2)
x_ax_postFR=[0:0.01:2];
hold on
    plot(x_ax_postFR,(rhoE.*x_ax_postFR.*(x_ax_postFR-cE))./tau_wEE,'Color',map(2,:),'LineWidth',width_of_lines);
    plot([cE,cE],[-1,1],':k','LineWidth',width_of_lines)
    plot(x_ax_postFR,zeros(length(x_ax_postFR),1),'black','LineWidth',width_of_lines)

    hold off
ylim([-3,4]*10^(-3))
XLABEL=xlabel('Postsyn. rate \nu^E in [Hz]');
YLABEL=ylabel('\Delta w^{EE}');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)


subplot(3,3,3)
hold on
for bb3=1:length(vec_nuI)
    plot(vec_rhoE,save_dwEE(bb3,:),'Color',map4(bb3,:),'LineWidth',width_of_lines)
end
plot(vec_rhoE,zeros(length(vec_rhoE),1),'black','LineWidth',width_of_lines)
hold off
xlim([0 3])
ylim([-10^(-3),2*10^(-3)])
XLABEL=xlabel('Presyn. exc. rate \rho^E in [Hz]');
YLABEL=ylabel('\Delta w^{EE}');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
hLegend=legend({},'FontSize',8,'FontName','Arial','location','northwest');
hLegend=legend('0','1','2','3','4','5','6','location', 'North' );
hLegend.ItemTokenSize = [3,18];
legend('boxoff')
set(gca,'linewidth',width_of_lines)

subplot(3,3,4)
x_ax_postFR=[0:0.01:4];
hold on
    plot(x_ax_postFR,(rhoE.*x_ax_postFR.*(x_ax_postFR-cE))./tau_wEE,'Color',map(2,:),'LineWidth',width_of_lines);
    plot(x_ax_postFR,((rhoI+wIE*rhoE).*(x_ax_postFR-cI))./tau_wEI,'Color',map2(2,:),'LineWidth',width_of_lines);
    plot(x_ax_postFR,zeros(length(x_ax_postFR),1),'black','LineWidth',width_of_lines)
    plot([cE,cE],[-1,1],':k','LineWidth',width_of_lines)
    hold off
ylim([-3,25]*10^(-3))
XLABEL=xlabel('Postsyn. rate \nu^E in [Hz]');
YLABEL=ylabel('\Delta w');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
hLegend=legend({'\Delta w^{EE}','\Delta w^{EI}'},'FontSize',8,'FontName','Arial','location','northwest');
hLegend.ItemTokenSize = [15,18];
legend('boxoff')
set(gca,'linewidth',width_of_lines)

subplot(3,3,5)
f=@(t,W) [(rhoE*max(NE*rhoE*W(1)-(rhoI+NI*wIE*rhoE)*W(2),0)*(max(NE*rhoE*W(1)-(rhoI+NI*wIE*rhoE)*W(2),0)-cE))/tau_wEE;
    ((rhoI+NI*wIE*rhoE)*(max(NE*rhoE*W(1)-(rhoI+NI*wIE*rhoE)*W(2),0)-cI))/tau_wEI];
wEI_0=[0.5,1,1.8];
wEE_0=[1.5,2.5,1.5];

size_grid=5;
w1 = linspace(0,size_grid,20);
w2 = linspace(0,size_grid,20);
[x,y] = meshgrid(w1,w2);

u = zeros(size(x));
v = zeros(size(x));
t=0; 
for i = 1:numel(x)
    Wprime = f(t,[x(i); y(i)]);
    u(i) = Wprime(1);
    v(i) = Wprime(2);
end

an=min(min(min(u(u<0))),min(min(v(v<0)))); % get the largest negative value
bn=max(max(max(u(u<0))),max(max(v(v<0)))); % get the smalles negative value
ap=min(min(min(u(u>0))),min(min(v(v>0)))); % get the smalles positive value
bp=max(max(max(u(u>0))),max(max(v(v>0)))); % get the largest positive value

cn=-0.2;
dn=-0.01;
cp=0.15; % chosen minimum pos value
dp=0.25; % chosen maximum pos value

v_old=v;

for i = 1:numel(x) %%%%%%%%%%%%%%%%%%%%%% this scales the length of the arrows so I can see them all
    if v(i)<0
        v(i) = cn+(dn-cn)/(bn-an)*(v(i)-an);
    elseif v(i)>0
        v(i) = cp+(dp-cp)/(bp-ap)*(v(i)-ap);        
    end
    mult_scale=v(i)/v_old(i);
        u(i) = u(i)*mult_scale;
end

hold on
quiver(x,y,u,v,'AutoScale','off','color',[0.5,0.5,0.5]);
xlim([0,max(w1)])
ylim([0,max(w2)])
for hh1=1:length(wEE_0)
    [ts,ys] = ode45(f,[0,50000],[wEE_0(hh1);wEI_0(hh1)]);
    plot(ys(:,1),ys(:,2),'Color',map3(6,:),'Linewidth',width_of_lines)
    scatter(wEE_0(hh1),wEI_0(hh1),'k','filled')
end

plot([0,size_grid],(NE*rhoE.*[0,size_grid])./NI./(rhoI+wIE*rhoE)-((rhoI+wIE*rhoE)*tau_wEE)./(NE*rhoE^2*tau_wEI),'k--','Linewidth',width_of_lines)
plot([0,size_grid],(NE*rhoE.*[0,size_grid]-cE)./NI./(rhoI+wIE*rhoE),'k','Linewidth',width_of_lines)
hold off
XLABEL=xlabel('w^{EE}');
YLABEL=ylabel('w^{EI}');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)

subplot(3,3,6)
hold on
plot([dt:dt:total_time]./1000,save_stuff2(:,3),'--','Color',map(2,:),'LineWidth',width_of_lines)
plot([dt:dt:total_time]./1000,save_stuff2(:,4),'--','Color',map2(2,:),'LineWidth',width_of_lines)
plot([dt:dt:total_time]./1000,save_stuff(:,3),'Color',map(2,:),'LineWidth',width_of_lines)
plot([dt:dt:total_time]./1000,save_stuff(:,4),'Color',map2(2,:),'LineWidth',width_of_lines)
hold off
ylim([0,6])
xlim([0,1])
XLABEL=xlabel('Time in [s]');
YLABEL=ylabel('w');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
hLegend=legend('w^{EE} unst.','w^{EI} unst.','w^{EE} sta.','w^{EI} sta.','location', 'North' );
hLegend.ItemTokenSize = [15,18];
legend('boxoff')
set(gca,'linewidth',width_of_lines)

  
%print(h1,'Figure_1_V2','-dpdf','-bestfit')