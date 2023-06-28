clear all 
close all



%==========================================================================
% Initialisation
%==========================================================================

mu = [0 0]; % mu=[mux mup]
q_R=10;
p_R=10;
N=51;
beta=0.05;
gamma=0.05;
delta=1;
t_final=10
t_0=1
h_step=100;

q_i=linspace(-q_R,q_R,N );
p_i=linspace(-p_R,p_R,N);
[qmesh,pmesh]=meshgrid(q_i,p_i);

%==========================================================================
% Conditioning via modified shooting method
%==========================================================================
Norm_hm=GetNormscapeAv(t_final,h_step,q_i,p_i,beta,gamma,delta);
CD=NaN(N,N);


figure
imagesc(q_i,p_i,Norm_hm)
title('quantum density')
colorbar
colormap(viridis)
set(gca,'YDir','normal')
xlabel('q')
ylabel('p')
% caxis([0 1])
c = colorbar('eastoutside');
