%%%%%% The following code generates S1 Fig. from Miehl & Gjorgjieva 2022
%%%%%% PLoS CB. https://doi.org/10.1371/journal.pcbi.1010682

wEE=1.5; % initial E-to-E weight strength
wEI=0.5; % initial I-to-E weight strength
wIE=0.5; % initial E-to-I weight strength
cE=1; % E postsynaptic LTD/LTP threshold 
cI=1; % I postsynaptic LTD/LTP threshold 

NE=1; % Number of presynaptic E neurons
NI=1;

rhoE=2; % Presynaptic E rate in [Hz]
rhoI=0.5; % External E rate onto I neurons in [Hz]
FR_E=max(NE*rhoE*wEE-NI*rhoI*wEI,0); % E postsynaptic firing rate in [Hz]
FR_I=rhoI+wIE*rhoE; % I firing rate in [Hz]


tau_FR_E=10; % Time constant for E neuron rate dynamics in [ms]
tau_FR_I=10; % Time constant for I neuron rate dynamics in [ms]
tau_wEE=1000; % Timescale for E plasticity in [ms]
tau_wEI=200; % Timescale for I plasticity in [ms]

map = brewermap(3,'Blues');
map2 = brewermap(3,'Reds');
map3 = brewermap(6,'Greens');

width_of_lines=1;
size_font=8;

h1=figure;
x_ax_preFR=[0:0.01:2];
y_ax_postFR=[0:0.01:2];

weight_change_matrix=(x_ax_preFR.*y_ax_postFR'.*(y_ax_postFR'-cE))./tau_wEE;
max_wcm=max(max(weight_change_matrix));
min_wcm=min(min(weight_change_matrix));

weight_change_matrix(weight_change_matrix<0)=weight_change_matrix(weight_change_matrix<0)./(-min_wcm);
weight_change_matrix(weight_change_matrix>0)=weight_change_matrix(weight_change_matrix>0)./(max_wcm);


vec=100:-20:0;
NNN = 128;
hex=['#ffffff','#d0d1e6','#a6bddb','#74a9cf','#2b8cbe','#045a8d']'; % Blue
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
map = interp1(vec,raw,linspace(100,0,NNN),'pchip');
hex2=['#ffffff','#fdd0a2','#fdae6b','#fd8d3c','#e6550d','#a63603']'; % orange
raw2 = sscanf(hex2','#%2x%2x%2x',[3,size(hex2,1)]).' / 255;
map2=interp1(vec,raw2,linspace(100,0,NNN),'pchip');
map3 = [flipud(map);map2];
colormap(map3)
imagesc(flipud(weight_change_matrix))

colorbar

xticks(0:50:200)
xticklabels(0:0.5:2)
yticks(0:50:200)
yticklabels(2:-0.5:0)

XLABEL=xlabel('Presyn. rate \rho^E');
YLABEL=ylabel('Postsyn. rate \nu^E');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)