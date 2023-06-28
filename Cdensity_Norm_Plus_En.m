clear all 
close all

%==========================================================================
% Initialisation
%==========================================================================

mu = [0 0]; % mu=[mux mup]
q_R=10;
p_R=10;
N=101;
beta=0.05;
gamma=0.05;
delta=1;
t_final=50;
h_step=0;
sigma_qp=pi
q_i=linspace(-q_R,q_R,N );
p_i=linspace(-p_R,p_R,N);
[qmesh,pmesh]=meshgrid(q_i,p_i);
dqs=abs(q_i(N)-q_i(N-1));
dps=abs(p_i(N)-p_i(N-1));
dqs=dqs*(1-1/N);
dps=dps*(1-1/N);


delta_upper_old=0.1;
delta_lower_old=0;
eps_stop=0.01; % small parameter %% h/2
itt_max=50;

%==========================================================================
% Load the quantum data
%==========================================================================

% parent_d = cd;    
% cd './FHus'
% Quant_dat = matfile('test_nefn_1_D51.mat');
% Quant_dat = matfile('cho_nefn2_D51');
% Hus=Quant_dat.Hus_av;
% Lh=0.5*(length(Hus)-1);
% nfq=sum(sum(Hus/Lh));
% cd(parent_d)
nefn=4
nfq=nefn*2*pi
%==========================================================================
% Conditioning via modified shooting method
%==========================================================================
Norm_hm=GetNormscapeCHO(t_final,h_step,q_i,p_i,beta,gamma,delta);
CD=NaN(N,N);


% 
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

% 
% for ittnormq
%     for ittnormp

for itt=1:itt_max
    % Reset arrays
    dnew=abs(delta_upper_old-delta_lower_old)/2; 
    delta_upper_new=delta_upper_old-dnew; % Move the right hand side 
    % Partition norm landscape
    eps_fwd=delta_upper_new;
    eps_bwd=eps_fwd;
    [Normp] = partition_loss(Norm_hm,eps_fwd,eps_bwd,'G');
    

    %======================================================================
    % The 2D smoothing
    %======================================================================
    dz=dqs*dps;
    CD=imgaussfilt(Normp,sigma_qp);

    nfc=sum(sum(CD*dqs*dps));
%      nfc=sum(sum(CD))
    
        figure(1)
        imagesc(q_i,p_i,CD)
        colorbar
        title('classical density')
        colormap(viridis)
        set(gca,'YDir','normal')
        xlabel('q')
        ylabel('p')
        caxis([0 1])
        c = colorbar('eastoutside');
%         return

        if abs(nfq-nfc)<eps_stop % Passes
                display('FINISHED ')
                break
        else % Fails and update interval
            if nfc>nfq
                delta_lower_old=delta_upper_new;
            elseif nfc<nfq
                delta_upper_old=delta_upper_new;
            end
        end

end

   
    %======================================================================
    %======================================================================

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



figure
imagesc(q_i,p_i,CD)
title('classical density')
colorbar
colormap(viridis)
set(gca,'YDir','normal')
xlabel('q')
ylabel('p')
caxis([0 1])
c = colorbar('eastoutside');


delta_upper_old=1;
delta_lower_old=0;
eps_stop=0.01; % small parameter %% h/2
itt_max=50;

