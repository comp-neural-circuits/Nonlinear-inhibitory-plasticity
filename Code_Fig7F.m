%%%%%% The following code generates Fig. 7F from Miehl & Gjorgjieva 2022
%%%%%% PLoS CB. https://doi.org/10.1371/journal.pcbi.1010682

close all
clear all

dt=0.1;
wEE_FF=0.06;
wEE_RC=0.001;
wEI_unsp=0.01;
wEI_spec=0.01;
wIE_unsp_FF=0.002;
wIE_spec_FF=0.001;
wIE_unsp_FB=0;
wIE_spec_FB=0;
cE=1;
cI=1;
rhoE=1;
rhoI_unsp=0;
rhoI_spec=0;

NE_RC=30;
NE_FF=40;
NI_unsp=20;
NI_spec=20;

pat_width_E=4;
stepsize_input=4; 

baseline_FR_E=rhoE*1;
pattern_FR_E=baseline_FR_E*4;

wEE_vec=ones(NE_FF,NE_RC).*wEE_FF.*rand(NE_FF,NE_RC).*3;
wEE_vec(wEE_vec<0)=0;
wEE_RC_mat=ones(NE_RC,NE_RC).*wEE_RC;
wEE_RC_mat=wEE_RC_mat-diag(diag(wEE_RC_mat)); 
wEI_unsp_vec=ones(NI_unsp,NE_RC).*wEI_unsp;
wEI_spec_vec=ones(NI_spec,NE_RC).*wEI_spec;
wIE_unsp_FF_mat=ones(NI_unsp,NE_FF*NE_RC).*wIE_unsp_FF; 

dummy_vec=zeros(1,NE_FF);
dummy_vec(1:pat_width_E)=ones(1,pat_width_E).*wIE_spec_FF;
wIE_spec_FF_mat_dummy=dummy_vec;
wIE_spec_FF_mat=[];

for uu=1:NI_spec-1
    if mod(uu,NI_spec/(NE_FF/pat_width_E))==0
        dummy_vec=circshift(dummy_vec,pat_width_E);
        wIE_spec_FF_mat_dummy=[wIE_spec_FF_mat_dummy;dummy_vec];
    else
       wIE_spec_FF_mat_dummy=[wIE_spec_FF_mat_dummy;dummy_vec];
    end
end

for uu2=1:NE_RC
    wIE_spec_FF_mat=[wIE_spec_FF_mat,wIE_spec_FF_mat_dummy];
end

wIE_unsp_FB_vec=ones(NI_unsp,NE_RC).*wIE_unsp_FB;
wIE_spec_FB_vec=ones(NI_spec,NE_RC).*wIE_spec_FB;
rhoE_vec=ones(NE_FF,NE_RC).*rhoE;
rhoE_vec2=ones(NE_FF*NE_RC,1)*rhoE;
rhoI_unsp_vec=ones(NI_unsp,1).*rhoI_unsp;
rhoI_spec_vec=ones(NI_spec,NE_RC).*rhoI_spec;

tau_FR_E=dt;
tau_FR_I=dt;
tau_wEE=1000;
tau_wEI=200;

total_time=2*200*1000;
time_pattern_pres=100;
pattern_start=1000;
pattern_end=2*200*1000;


start_disinh2=30*1000;
end_disinh2=200*1000;

save_time=100;

FR_E=max(diag(rhoE_vec'*wEE_vec)-(rhoI_unsp_vec'*wEI_unsp_vec)'-diag(rhoI_spec_vec'*wEI_spec_vec),0);
FR_I_unsp_vec=max(rhoI_unsp_vec+wIE_unsp_FF_mat*rhoE_vec2+wIE_unsp_FB_vec*FR_E,0);
FR_I_spec_vec=max(rhoI_spec_vec+wIE_spec_FF_mat*rhoE_vec2+wIE_spec_FB_vec*FR_E,0);
counter=0;

counter_pattern=0;

circular_vec=[1:NE_FF,1:NE_FF,1:NE_FF];

choose_E=zeros(1,pat_width_E);
vec_of_pattern_E=zeros(1,pat_width_E);

slow_recovery=linspace(-2,0,100000);
counter_recov=0;

for tt=dt:dt:total_time
    
    if mod(round(tt,2),time_pattern_pres)==0 && tt>pattern_start && tt<pattern_end
        xx=randi(NE_FF/stepsize_input)*stepsize_input;
        
               
        input_FR_E_vec=baseline_FR_E.*ones(1,NE_FF);
        choose_E(1,:)=NE_FF+xx+1:NE_FF+xx+pat_width_E;
        vec_of_pattern_E(1,:)=circular_vec(choose_E(1,:));
        dummy_input=input_FR_E_vec(1,:);
        

        dummy_input(vec_of_pattern_E(1,:))=pattern_FR_E;

        input_FR_E_vec(1,:)=dummy_input;
        
        dummy_input_mat=[];
        dummy_input_mat2=[];
        for uu3=1:NE_RC
           dummy_input_mat=[dummy_input_mat,input_FR_E_vec']; 
           dummy_input_mat2=[dummy_input_mat2;input_FR_E_vec'];
        end
        rhoE_vec=dummy_input_mat;
        rhoE_vec2=dummy_input_mat2;
        
        counter_pattern=counter_pattern+1;
        
    end

    
    if round(tt,2)==start_disinh2
        rhoI_spec_vec=ones(NI_spec,NE_RC).*(-2);
    elseif round(tt,2)>=end_disinh2 && counter_recov<length(slow_recovery)
        
        counter_recov=counter_recov+1;
        rhoI_spec_vec=ones(NI_spec,NE_RC).*slow_recovery(counter_recov);
    end

    
    FR_I_unsp_vec=FR_I_unsp_vec+(-FR_I_unsp_vec+max(rhoI_unsp_vec+wIE_unsp_FF_mat*rhoE_vec2+wIE_unsp_FB_vec*FR_E,0))./tau_FR_I*dt;
    FR_I_spec_vec=FR_I_spec_vec+(-FR_I_spec_vec+max(rhoI_spec_vec+wIE_spec_FF_mat*rhoE_vec2+wIE_spec_FB_vec*FR_E,0))./tau_FR_I*dt;
    FR_E=FR_E+(-FR_E+max(diag(rhoE_vec'*wEE_vec)+wEE_RC_mat*FR_E-(FR_I_unsp_vec'*wEI_unsp_vec)'-diag(FR_I_spec_vec'*wEI_spec_vec),0))./tau_FR_E*dt;

    wEE_vec=wEE_vec+(rhoE_vec.*repmat((FR_E.*(FR_E-cE))',NE_FF,1))./tau_wEE*dt;
    wEE_RC_mat=wEE_RC_mat+(FR_E*(FR_E.*(FR_E-cE))')./tau_wEE*dt;
    wEE_RC_mat=wEE_RC_mat-diag(diag(wEE_RC_mat));
    
    wEI_unsp_vec=wEI_unsp_vec+(FR_I_unsp_vec*(FR_E.*(FR_E-cI))')./tau_wEI*dt;
    wEI_spec_vec=wEI_spec_vec+(FR_I_spec_vec.*repmat((FR_E.*(FR_E-cI))',NI_unsp,1))./tau_wEI*dt;
    wEE_vec(wEE_vec<0)=0;
    wEE_RC_mat(wEE_RC_mat<0)=0;
    wEI_unsp_vec(wEI_unsp_vec<0)=0;
    wEI_spec_vec(wEI_spec_vec<0)=0;

    if mod(round(tt,2),save_time)==0
        counter=counter+1;
        save_stuff_FRE(counter,:)=FR_E;
        save_stuff_FRI_unsp(counter,1)=FR_I_unsp_vec(1,1);
        save_stuff_FRI_spec(counter,1)=FR_I_spec_vec(1,1);
        
        for oo=1:NE_RC
            save_stuff_wEE(counter,1+NE_FF*(oo-1):NE_FF*oo)=wEE_vec(:,oo);
            save_stuff_wEI_unsp(counter,1+NI_unsp*(oo-1):NI_unsp*oo)=wEI_unsp_vec(:,oo);
            save_stuff_wEI_spec(counter,1+NI_spec*(oo-1):NI_spec*oo)=wEI_spec_vec(:,oo);
        end
        
        save_stuff_wEE_RC(counter,1)=wEE_RC_mat(2,1);
        save_stuff_wEE_RC(counter,2)=wEE_RC_mat(1,2);
        save_stuff_wEE_RC(counter,3)=wEE_RC_mat(3,2);
        save_stuff_wEE_RC(counter,4)=wEE_RC_mat(2,3);
        save_stuff_wEE_RC(counter,5)=wEE_RC_mat(1,3);
        save_stuff_wEE_RC(counter,6)=wEE_RC_mat(3,1);
            
    end
    
end


%% Plot figures
map2 = brewermap(3,'Set1');

width_of_lines=1;
size_font=8;

h2=figure;

%%%%%% make network plot here
subplot(3,3,[5 6 7 8])
s=[];
t=[];
w_in_vec=[];

wEE_RC_mat_cut=wEE_RC_mat;
wEE_RC_mat_cut(wEE_RC_mat_cut<0.03)=0;

for hh=1:NE_RC
    xx3=[];
    s=[s,ones(1,NE_RC-1)*hh];
    if hh>1
        xx3=[xx3,1:hh-1];
    end
    xx3=[xx3,hh+1:NE_RC];
t=[t,xx3];
w_in_vec=[w_in_vec,wEE_RC_mat_cut(xx3,hh)'];
end
w_in_vec=round(w_in_vec,2);

counter=0;
for uu=1:length(w_in_vec)
    if w_in_vec(uu)~=0
        counter=counter+1;
        s2(counter)=s(uu);
        t2(counter)=t(uu);
        w_in_vec2(counter)=w_in_vec(uu);
    end
end

for nn2=1:NE_RC
    save_vals=[];
    for nn3=1:round(NE_FF/stepsize_input)
        save_vals(nn3)=mean(wEE_vec(1+stepsize_input*(nn3-1):stepsize_input*nn3,nn2));
    end

[a,b]=max(save_vals);
save_max(nn2,1:2)=[a,b];

end


G = digraph(s2,t2,w_in_vec2);
LWidths = 2*G.Edges.Weight/max(G.Edges.Weight);
p=plot(G,'LineWidth',LWidths);
p=plot(G,'Nodelabel',save_max(:,2),'LineWidth',LWidths);
%p=plot(G,'Nodelabel',save_max(:,2),'Edgelabel',G.Edges.Weight,'LineWidth',LWidths);
p.Marker = '^';
p.NodeColor = 'b';
p.EdgeColor = '[0.7,0.7,0.7]';
set(gca,'XColor', 'none','YColor','none')

h2.Renderer='Painters';
%print(h2,'Fig6E_1','-dpdf','-bestfit')

%%
h4=figure;
map8=brewermap(10,'Greens');

[sorted,indx]=sort(save_max(:,2),'ascend');
wEE_RC_row_sorted=wEE_RC_mat(indx,:);
wEE_RC_row_sorted=wEE_RC_row_sorted(:,indx);

subplot(3,3,7)
imagesc(sorted)
colormap(map8)
colorbar

h3=figure;

subplot(3,3,8)
imagesc(wEE_vec)
colorbar
bla=gray;
bla=flipud(bla);
colormap(bla)
axis square
XLABEL=xlabel('RC exc. neuron');
YLABEL=ylabel('FF exc. neuron');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)


subplot(3,3,9)
imagesc(wEE_RC_row_sorted)
colorbar
bla=gray;
bla=flipud(bla);
colormap(bla)
axis square
XLABEL=xlabel('RC exc. neuron - post');
YLABEL=ylabel('RC exc. neuron - pre');
set([XLABEL,YLABEL],'FontName','Arial');
set([XLABEL,XLABEL],'FontSize',size_font);
set(gca,'FontSize',size_font,'FontName','Arial');
set(gca,'linewidth',width_of_lines)

h3.Renderer='Painters';
%print(h3,'Fig6E_2','-dpdf','-bestfit')