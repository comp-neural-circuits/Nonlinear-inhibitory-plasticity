%%%%%% The following code generates Fig. 5 from Miehl and Gjorgjieva 2021
close all
clear all

%% Parameter definitions
nr_vals=1000;
wEI_vec=3*rand(nr_vals,1);
wEE_vec=3*rand(nr_vals,1);

dt=0.1;
wIE=0.5;
cE=1;
cI=1;
rhoI=0.5;

% these are the examples which I color in the "data" panals
wEI_0=[0.5,1,1.8];
wEE_0=[1.5,2.5,1.5];

wEI_vec=[wEI_vec;wEI_0'];
wEE_vec=[wEE_vec;wEE_0'];

max_EIbalance=12;

NE=1;
NI=1;

tau_FR_E=10;
tau_FR_I=10;
tau_wEE=1000;
tau_wEI=200;

total_time=100;

rhoE_0=2;
rhoE=rhoE_0;

counter=0;

%% Simulation start
for xx=1:nr_vals+length(wEI_0)
    wEE=wEE_vec(xx);
    wEI=wEI_vec(xx);
    
    FR_E=max(NE*rhoE*wEE-NI*rhoI*wEI,0);
    FR_I=rhoI+wIE*rhoE;
    
    counter=counter+1;
    
    save_stuff(counter,1)=wEE;
    save_stuff(counter,2)=wEI;
    
    if wEE/wEI<max_EIbalance

    for tt=dt:dt:total_time


        FR_E=FR_E+(-FR_E+max(NE*rhoE*wEE-NI*FR_I*wEI,0))/tau_FR_E*dt;
        FR_I=FR_I+(-FR_I+rhoI+wIE*rhoE)/tau_FR_I*dt;

        wEE=wEE+(rhoE*FR_E*(FR_E-cE))/tau_wEE*dt;
        wEI=wEI+(FR_I*FR_E*(FR_E-cI))/tau_wEI*dt;
        wEE(wEE<0)=0;
        wEI(wEI<0)=0;
    end

    if wEE/wEI<max_EIbalance
        save_stuff(counter,3)=wEE;
        save_stuff(counter,4)=wEI;
    else
        save_stuff(counter,1)=NaN;
        save_stuff(counter,2)=NaN;
        save_stuff(counter,3)=NaN;
        save_stuff(counter,4)=NaN;        
    end
        
        
    else
        save_stuff(counter,1)=NaN;
        save_stuff(counter,2)=NaN;
        save_stuff(counter,3)=NaN;
        save_stuff(counter,4)=NaN;
    end
end

%% Plot figures
upperlim=2/tau_wEE;
lowerlim=-1.1/tau_wEE;

map = brewermap(3,'Blues');
map0 = brewermap(4,'Set1');
map2 = brewermap(3,'Reds');
map3 = brewermap(6,'Greens');

width_of_lines=1;
size_font=8;

h2=figure;

xx=0.1:0.01:5;  
EIbal=(NI.*(rhoI+wIE*rhoE_0).*xx+cE)./NE./(rhoE_0.*xx);

subplot(3,3,7)
hold on
plot(xx,EIbal,'k','LineWidth',width_of_lines)
plot(xx,ones(length(xx),1).*NI.*(rhoI+wIE*rhoE_0)./NE./rhoE_0,'--k','LineWidth',width_of_lines)
hold off
xlim([0 max(xx)])
ylim([0 5])
XLABEL=xlabel('w^{EI}');
YLABEL=ylabel('R^{E/I}');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)


subplot(3,3,1)
rhoE_vec=0.01:0.01:5;
EIbal2=(rhoI+wIE*rhoE_vec)./(rhoE_vec);
hold on
plot(rhoE_vec,EIbal2,'k','LineWidth',width_of_lines)
hold off
xlim([0 max(rhoE_vec)])
ylim([0 6])
XLABEL=xlabel('\rho^{E}');
YLABEL=ylabel('R^{E/I}');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)





% phase plane plot here
subplot(3,3,2)
f=@(t,W) [(rhoE_0*max(NE*rhoE_0*W(1)-NI*(rhoI+wIE*rhoE_0)*W(2),0)*(max(NE*rhoE_0*W(1)-NI*(rhoI+wIE*rhoE_0)*W(2),0)-cE))/tau_wEE;
    ((rhoI+wIE*rhoE_0)*max(NE*rhoE_0*W(1)-NI*(rhoI+wIE*rhoE_0)*W(2),0)*(max(NE*rhoE_0*W(1)-NI*(rhoI+wIE*rhoE_0)*W(2),0)-cI))/tau_wEI];


size_grid=5;
w1 = linspace(0,size_grid,20);
w2 = linspace(0,size_grid,20);
[x,y] = meshgrid(w1,w2);

u = zeros(size(x));
v = zeros(size(x));
t=0; % we want the derivatives at each point at t=0, i.e. the starting time
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
cp=0.05; % chosen minimum pos value
dp=0.3; % chosen maximum pos value

u_old=u;

for i = 1:numel(x)

    if u(i)<0
        u(i) = cn+(dn-cn)/(bn-an)*(u(i)-an);
    elseif u(i)>0
        u(i) = cp+(dp-cp)/(bp-ap)*(u(i)-ap);        
    end
    mult_scale=u(i)/u_old(i);
        v(i) = v(i)*mult_scale;

end


hold on
quiver(x,y,u,v,'AutoScale','off','color',[0.5,0.5,0.5]);
xlim([0,max(w1)])
ylim([0,max(w2)])
symbols=['x','x','x'];
for hh1=1:length(wEE_0)
    [ts,ys] = ode45(f,[0,100],[wEE_0(hh1);wEI_0(hh1)]);
    plot(ys(:,1),ys(:,2),'Color',map3(6,:),'Linewidth',width_of_lines)
    scatter(wEE_0(hh1),wEI_0(hh1),'k','filled')
    scatter(ys(end,1),ys(end,2),'k',symbols(hh1))
end

plot([0,size_grid],(NE*rhoE_0.*[0,size_grid])./NI./(rhoI+wIE*rhoE_0),'color',[0.7,0.7,0.7],'Linewidth',width_of_lines)
plot([0,size_grid],(NE*rhoE_0.*[0,size_grid]-cE)./NI./(rhoI+wIE*rhoE_0),'k','Linewidth',width_of_lines)
hold off
XLABEL=xlabel('w^{EE}');
YLABEL=ylabel('w^{EI}');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)



subplot(3,3,3)
hold on
scatter(save_stuff(:,1)./save_stuff(:,2),save_stuff(:,3)./save_stuff(:,4),10,'k','filled')
scatter(save_stuff(end-length(wEI_0)+1:end,1)./save_stuff(end-length(wEI_0)+1:end,2),save_stuff(end-length(wEI_0)+1:end,3)./save_stuff(end-length(wEI_0)+1:end,4),30,'r','x')
plot([0,max([save_stuff(:,1)./save_stuff(:,2);save_stuff(:,3)./save_stuff(:,3)])],[0,max([save_stuff(:,1)./save_stuff(:,2);save_stuff(:,3)./save_stuff(:,3)])],'--','color',[0.6,0.6,0.6],'LineWidth',width_of_lines)
plot([0,max_EIbalance],[(rhoI+wIE*rhoE)/rhoE,(rhoI+wIE*rhoE)/rhoE],'color',[0.6,0.6,0.6],'LineWidth',width_of_lines)
hold off
xlim([0,max_EIbalance])
XLABEL=xlabel('E/I ratio before');
YLABEL=ylabel('E/I ratio after');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)


subplot(3,3,4)
hold on
scatter(100./save_stuff(:,1).*save_stuff(:,3),100./save_stuff(:,2).*save_stuff(:,4),10,'k','filled')
scatter(100./save_stuff(end-length(wEI_0)+1:end,1).*save_stuff(end-length(wEI_0)+1:end,3),100./save_stuff(end-length(wEI_0)+1:end,2).*save_stuff(end-length(wEI_0)+1:end,4),30,'r','x')
plot([0,max(100./save_stuff(:,1).*save_stuff(:,3))],[100,100],'--','color',[0.6,0.6,0.6],'LineWidth',width_of_lines)
hold off
xlim([80,max(100./save_stuff(:,1).*save_stuff(:,3))])
XLABEL=xlabel('\Delta Excitation (%)');
YLABEL=ylabel('\Delta Inhibition (%)');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)



subplot(3,3,5)
hold on
scatter(save_stuff(:,1)./save_stuff(:,2),100./save_stuff(:,1).*save_stuff(:,3),10,'k','filled')
scatter(save_stuff(end-length(wEI_0)+1:end,1)./save_stuff(end-length(wEI_0)+1:end,2),100./save_stuff(end-length(wEI_0)+1:end,1).*save_stuff(end-length(wEI_0)+1:end,3),30,'r','x')
plot([0,max_EIbalance],[100,100],'--','color',[0.6,0.6,0.6],'LineWidth',width_of_lines)
hold off
xlim([0,max_EIbalance])
XLABEL=xlabel('E/I ratio before');
YLABEL=ylabel('\Delta Excitation (%)');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)



subplot(3,3,6)
hold on
scatter(save_stuff(:,1)./save_stuff(:,2),100./save_stuff(:,2).*save_stuff(:,4),10,'k','filled')
scatter(save_stuff(end-length(wEI_0)+1:end,1)./save_stuff(end-length(wEI_0)+1:end,2),100./save_stuff(end-length(wEI_0)+1:end,2).*save_stuff(end-length(wEI_0)+1:end,4),30,'r','x')
plot([0,max_EIbalance],[100,100],'--','color',[0.6,0.6,0.6],'LineWidth',width_of_lines)
hold off
xlim([0,max_EIbalance])
XLABEL=xlabel('E/I ratio before');
YLABEL=ylabel('\Delta Inhibition (%)');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)

h2.Renderer='Painters';
%print(h2,'Fig_5','-dpdf','-bestfit')