%%%%%% The following code generates S2B,C Fig. from Miehl & Gjorgjieva 2022
%%%%%% PLoS CB. https://doi.org/10.1371/journal.pcbi.1010682

clear all
vec_rhoE=[0:0.01:5];

dt=0.1;
wEI=0.5;
wIE_FF=1;
wIE_FB=1;
cE=1;
cI=1;
wEE=1.5;
rhoI=0.5;
wEI2=wEI;
wIE2_FF=0.4;
wIE2_FB=0.4;
wEE2=wEE;
rhoI2=rhoI;

NE=1;

tau_FR_E=10;
tau_FR_I=10;
tau_wEE=500;
tau_wEI=1000;

total_time=3000;

rhoE_0=2;
rhoE=rhoE_0;
rhoE2=rhoE_0;
rhoE3=rhoE_0;
rhoE4=rhoE_0;

FR_E=max(rhoE*wEE-rhoI*wEI,0);
FR_I=rhoI+wIE_FF*rhoE+wIE_FB*FR_E;
FR_E2=max(rhoE2*wEE2-rhoI2*wEI2,0);
FR_I2=rhoI2+wIE2_FF*rhoE2+wIE2_FB*FR_E2;
counter=0;

for tt=dt:dt:total_time
    
    counter=counter+1;

    FR_E=FR_E+(-FR_E+max(rhoE*wEE-FR_I*wEI,0))/tau_FR_E*dt;
    FR_I=FR_I+(-FR_I+rhoI+wIE_FF*rhoE+wIE_FB*FR_E)/tau_FR_I*dt;

    wEE=wEE+(rhoE*FR_E*(FR_E-cE))/tau_wEE*dt;
    wEI=wEI+(FR_I*FR_E*(FR_E-cI))/tau_wEI*dt;
    wEE(wEE<0)=0;
    wEI(wEI<0)=0;

    save_stuff(counter,1)=FR_E;
    save_stuff(counter,2)=FR_I;
    save_stuff(counter,3)=wEE;
    save_stuff(counter,4)=wEI;
    
end


map = brewermap(3,'Blues');
map2 = brewermap(3,'Reds');


width_of_lines=1.5;
size_font=8;

h1=figure;

% this is the plot of the plasticity rules
subplot(3,3,2)
x_ax_postFR=[0:0.01:2];
hold on
    plot(x_ax_postFR,(rhoE_0.*x_ax_postFR.*(x_ax_postFR-cE))./tau_wEE,'Color',map(2,:),'LineWidth',width_of_lines);
    plot(x_ax_postFR,((rhoI+wIE_FF*rhoE_0+wIE_FB*FR_E).*x_ax_postFR.*(x_ax_postFR-cI))./tau_wEI,'Color',map2(2,:),'LineWidth',width_of_lines);
    plot(x_ax_postFR,zeros(length(x_ax_postFR),1),'black','LineWidth',width_of_lines)
    plot([cE,cE],[-1,1],':k','LineWidth',width_of_lines)
    hold off
ylim([-3,4]*10^(-3))
XLABEL=xlabel('Postsyn. rate \nu^E in [Hz]');
YLABEL=ylabel('\Delta w');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
hLegend=legend({'\Delta w^{EE}','\Delta w^{EI}'},'FontSize',8,'FontName','Arial','location','northwest');
hLegend.ItemTokenSize = [15,18];
legend('boxoff')
set(gca,'linewidth',width_of_lines)



subplot(3,3,4)
yyaxis right
hold on
plot([dt:dt:total_time]./1000,save_stuff(:,1));
plot([dt:dt:total_time]./1000,save_stuff(:,2));
hold off
ylim([0,6])
XLABEL=xlabel('Time in [s]');
YLABEL=ylabel('Firing rate in [Hz]');


yyaxis left
hold on
plot([dt:dt:total_time]./1000,save_stuff(:,3));
plot([dt:dt:total_time]./1000,save_stuff(:,4));
hold off
ylim([0,1.7])
XLABEL=xlabel('Time in [s]');
YLABEL=ylabel('w');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)

