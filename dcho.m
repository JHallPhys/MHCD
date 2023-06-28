function pdr=dcho(t,y,beta,gamma,delta) 



q=y(1);
p=y(2);


pdr=[p-beta*p^3-beta*p*q^2-gamma*q;-q+beta*q^3+beta*q*p^2+delta-gamma*p;(-2*gamma*(q^2+p^2))*y(3)];
% pdr=[y(2)-gamma*y(1);-y(1)-gamma*y(2);-gamma*(y(1)^2+y(2)^2)*y(3)];



% pdr=[y(2);-y(1);1];
% pdr=[p-gamma*q;-q-gamma*p];



end