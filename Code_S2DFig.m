%%%%%% The following code generates S2D Fig. from Miehl & Gjorgjieva 2022
%%%%%% PLoS CB. https://doi.org/10.1371/journal.pcbi.1010682


close all
clear all

vec_wIE_FF=[0:0.05:1.5];
vec_wIE_FB=[0:0.05:1.5];

save_stability=zeros(length(vec_wIE_FB),length(vec_wIE_FF));
save_stability_before=zeros(length(vec_wIE_FB),length(vec_wIE_FF));


for xx=1:length(vec_wIE_FF)
    for xx2=1:length(vec_wIE_FB)

dt=0.1;
wEI=0.5;
wIE_FB=vec_wIE_FB(xx2);
wIE_FF=vec_wIE_FF(xx);
cE=1;
cI=1;
wEE=1.5;
rhoI=0.5;


NE=1;

tau_FR_E=10;
tau_FR_I=10;
tau_wEE=500;
tau_wEI=1000;

total_time=10000;

rhoE=2; % FP I want to reach first


FR_E=max(rhoE*wEE-rhoI*wEI,0);
FR_I=rhoI+wIE_FB*FR_E+wIE_FF*rhoE;
counter=0;

for tt=dt:dt:total_time
    
    
    counter=counter+1;

    FR_E=FR_E+(-FR_E+max(rhoE*wEE-FR_I*wEI,0))/tau_FR_E*dt;
    FR_I=FR_I+(-FR_I+rhoI+wIE_FB*FR_E+rhoE*wIE_FF)/tau_FR_I*dt;

    wEE_old=wEE;
    wEE=wEE+(rhoE*FR_E*(FR_E-cE))/tau_wEE*dt;
    wEI=wEI+(FR_I*FR_E*(FR_E-cI))/tau_wEI*dt;
    wEE(wEE<0)=0;
    wEI(wEI<0)=0;
    

end
save_stability(xx2,xx)=wEE-wEE_old;
save_endweights(xx2,xx)=wEE;
    end
end
%%
stable_or_not=zeros(length(vec_wIE_FB),length(vec_wIE_FF));
for gg=1:length(vec_wIE_FF)
    for gg2=1:length(vec_wIE_FB)
        if save_stability(gg2,gg)>10^(-5) || isnan(save_stability(gg2,gg))==1
            stable_or_not(gg2,gg)=1;
        else
            stable_or_not(gg2,gg)=0;
        end
    end
end
%%
h3=figure; hAxes=gca;
imagesc(hAxes,stable_or_not)
colormap(hAxes,[0.5,0.5,0.5;1,1,1])
axis square
set(gca,'XTick',[0:8:32],'XTickLabel',[0:0.4:1.6])
set(gca,'YTick',[0:8:32],'YTickLabel',[0:0.4:1.6])
XLABEL=xlabel('w^{IE}_{FF}');
YLABEL=ylabel('w^{IE}_{FB}');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',16);
set(gca,'FontSize',16,'FontName','Arial');
set(h3,'Units','Inches');
set(gcf, 'PaperPositionMode','auto');
name2 = sprintf('stability_FF_vs_FB_V2.pdf');
print(h3,name2,'-dpdf','-r500');  
