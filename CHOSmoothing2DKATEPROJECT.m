clear all 
close all



%==========================================================================
% Initialisation
%==========================================================================

mu = [0 0]; % mu=[mux mup]
q_R=20;
p_R=20;
N=101;
beta=0.01;
gamma=0.01;
delta=1;
t_final=20;
sigma_qp=pi
q_i=linspace(-q_R,q_R,N );
p_i=linspace(-p_R,p_R,N);
[qmesh,pmesh]=meshgrid(q_i,p_i);
dqs=abs(q_i(N)-q_i(N-1));
dps=abs(p_i(N)-p_i(N-1));
dqs=dqs*(1-1/N);
dps=dps*(1-1/N);
nefn=30
h_step=20;

delta_upper_old=1;
delta_lower_old=0;
eps_stop=0.01; % small parameter %% h/2
itt_max=50;

%==========================================================================
% Load the quantum data
%==========================================================================

% parent_d = cd;    
% cd './FHus'
% Quant_dat = matfile('mhat_nefn_4.mat');
% % Quant_dat = matfile('testit2.mat');
% % Quant_dat = matfile('cho_nefn2_D51');
% Hus=Quant_dat.Hus_av;
% % Lh=0.5*(length(Hus)-1);
% nfq=sum(sum(Hus))/N;
% cd(parent_d)

nfq=nefn*2*pi
%==========================================================================
% Conditioning via modified shooting method
%==========================================================================
Norm_hm=GetNormscapeCHO(t_final,h_step,q_i,p_i,beta,gamma,delta);
CD=NaN(N,N);


figure
imagesc(q_i,p_i,Norm_hm)
colorbar
colormap(viridis)
set(gca,'YDir','normal')
xlabel('q')
ylabel('p')
% caxis([0 1])
c = colorbar('eastoutside');

% figure
% imagesc(q_i,p_i,log(Norm_hm))
% % surf(q_i,p_i,CD)
% % title('Norm')
% colorbar
% colormap(viridis)
% set(gca,'YDir','normal')
% xlabel('q')
% ylabel('p')
% % caxis([0 1])
% c = colorbar('eastoutside');
% 
% 
% return
% 
% for ittnormq
%     for ittnormp

for itt=1:itt_max
    itt
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

    nfc=sum(sum(CD*dqs*dps))
%      nfc=sum(sum(CD))
    
        figure(2)
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

    % 
% figure
% imagesc(q_i,p_i,Norm_hm)
% % surf(q_i,p_i,CD)
% title('Norm')
% colorbar
% colormap(viridis)
% set(gca,'YDir','normal')
% xlabel('q')
% ylabel('p')
% % caxis([0 1])
% c = colorbar('eastoutside');
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% CD=Norm_hm
% 
% figure
% imagesc(q_i,p_i,Hus)
% title('quantum density')
% colorbar
% colormap(viridis)
% set(gca,'YDir','normal')
% xlabel('q')
% ylabel('p')
% caxis([0 1])
% c = colorbar('eastoutside');

% 
% figure
% surf(q_i,p_i,CD)
% title('classical density')
% colorbar
% colormap(viridis)
% set(gca,'YDir','normal')
% xlabel('q')
% ylabel('p')
% caxis([0 1])
% c = colorbar('eastoutside');

% return

% [~,qind]=max(max(Hus));
% [~,cind]=max(max(CD));
% % [~,nind]=max(max(CD));
% figure
% hold on 
% plot(linspace(-q_R,q_R,length(Hus(:,qind))),Hus(:,qind),'b.-','Markersize',5)
% plot(q_i,CD(:,cind),'r.-','Markersize',5)
% legend('Quantum','Classical')

% 
% a=CD(:,cind)*1.2523
% b = normcdf(a,0,1/sqrt(2))
% 
% 
% figure
% hold on
% plot(linspace(-q_R,q_R,length(Hus(:,qind))),Hus(:,qind),'b.-','Markersize',5)
% plot(q_i,CD(:,cind),'r.-','Markersize',5)
% 
