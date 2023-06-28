function [Norm_out]=GetNormscapeAv(t_final,h_step,q_i,p_i,beta,gamma,delta)

%==========================================================================
% ODE Setup
%==========================================================================

t_i=0; % Start time
t_f=t_final; % Final time
% t_split=linspace(t_i,t_f,h_step); % Interval to integrate over 

z=zeros(3,1); % Vector of initial conditions z=[q(0),p(0),norm(q,p,t=0)]'

%==========================================================================
% Populate the Norm landscape
%==========================================================================


reverseStr = '';

for j=1:length(q_i)

    % Counter for the number we are on
    msg = sprintf('GetNorm %d/%d', j, length(q_i));
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));


    for k = 1:length(p_i)
        
        % 1) Initial conditions in phase space       
        z(1)=q_i(j);
        z(2)=p_i(k);
        z(3)=1;
        
          % 2) Propagate Norm
        
        [t,y] = ode89(@(t,y) dcho(t,y,beta,gamma,delta) ,linspace(0, t_final,h_step),z); % Integrate
        
        % 3) Save value of norm averaged over motion
     
        Norm_out(k,j)=mean(y(:,3));
%         Norm_out(k,j)=y(end,3);
    end
       
    
end    
    

%==========================================================================
% Populate the Norm landscape
%==========================================================================





end