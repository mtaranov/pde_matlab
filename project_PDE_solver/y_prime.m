function dydt = y_prime(t,y)
load n
ne=2;
dc=10;

% 
% konRtot=1.32;
% koff=10^-6;
% kdeg=2*10^-4;
% V_Rtot=5*10^-4;
% 
% konRtot=0.01;
% koff=10^-6;
% kdeg=2*10^-4;
% V_Rtot=5*10^-4;

% konRtot=1.32;
% koff=10^-6;
% kdeg=2*10^-4;
% V_Rtot=5*10^-5;
% 
konRtot=0.01;
koff=10^-6;
kdeg=2*10^-4;
V_Rtot=5*10^-5;

save('konRtot','konRtot')
save('koff','koff')
save('kdeg','kdeg')
save('V_Rtot','V_Rtot')

dydt=zeros(ne*n,1);
dh=100/n;%space step

k=1;
dydt(k)=V_Rtot-konRtot*y(k)*(1-y(k+1))+koff*y(k+1);
    
for k=3:2:ne*n-2
    L = k - ne;
    R = k + ne;
    dydt(k) = dc*(y(R)-2*y(k)+y(L))/dh^2 -konRtot*y(k)*(1-y(k+1))+koff*y(k+1);
end

dydt(k+2)= 0; 

for k2=2:2:ne*n-2
    dydt(k2)=konRtot*y(k2-1)*(1-y(k2))-koff*y(k2)-kdeg*y(k2);
end
k2=k2+2;

dydt(k2)=0;









