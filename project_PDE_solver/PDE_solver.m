% PDE solver using 2nd order RK

clean;
to=0;
n_curves=3;%number of successive curves
t_curves=600;%time btw succesive curves; A-300, B-600, C,D-1800
tmax=t_curves*n_curves;%final time; 
ne=2;

%ORDER OF ACCURACY IN TIME

n=20; %number of space steps
save('n','n')
ne=2;
co=zeros(1,ne*n);

    dt=0.8;
    m=2*tmax/dt;
    a=zeros(3,n*ne,m+1);
    y=zeros(ne*n,1);
    for i=1:3
       dt=dt/2;
       t=0:dt:tmax;
       n_time=tmax/dt; %number of time steps
  
            for j=1:n_time
            k1=dt*y_prime(t(j),y(:,j));
            k2=dt*y_prime(t(j)+dt/2,y(:,j)+k1/2);
            y(:,j+1)=y(:,j)+k2;
        end
        if i==1
            a(1,:,:)=y;
        elseif i==2
            a(2,:,:)=y(:,1:2:n_time+1);
        elseif i==3
            a(3,:,:)=y(:,1:4:n_time+1);
        end
    
    end


       p1=log2(abs((a(1,1:2:38,:)-a(2,1:2:38,:))./(a(2,1:2:38,:)-a(3,1:2:38,:))));
       p2=log2(abs((a(2,2:2:38,:)-a(2,2:2:38,:))./(a(2,2:2:38,:)-a(3,2:2:38,:))));


load V_Rtot
load konRtot
load koff
load kdeg

 for i=1:19
     for j=1:1501
         pnew1(i,j)=p1(1,i,j);
         
     end
 end
hold on; 

for i=1:19
     plot(pnew1(i,:),'r');
 end


xlabel('time')
ylabel('Order of accuracy, P1')
title([' nu/ R_t_o_t = ',num2str(V_Rtot),'   k_o_nR_t_o_t = ',num2str(konRtot), '  k_o_f_f = ',num2str(koff)])


    %ORDER OF ACCURACY IN SPACE

    n=40;%number of space steps
    dt=1;%time step
    t=0:dt:tmax;
    n_time=tmax/dt;%number of time steps
    a1=[];
    a2=[];
    a3=[];
     for i=1:3
                n=n/2;
                y=zeros(ne*n,1);
                save('n','n')
            for j=1:n_time
                k1=dt*y_prime(t(j),y(:,j));
                k2=dt*y_prime(t(j)+dt/2,y(:,j)+k1/2);
                y(:,j+1)=y(:,j)+k2;
             end
                 if i==1
                    a3=y;
                     elseif i==2
                       a2=y;
                     elseif i==3
                       a1=y;
                 end
     end
 
     
 p1=log2(abs((a1(3:2:10,:)-a2(5:4:20,:))./(a2(5:4:20,:)-a3(9:8:40,:))));
 p2=log2(abs((a1(4:2:10,:)-a2(6:4:20,:))./(a2(6:4:20,:)-a3(10:8:40,:))));

 
 plot(1:1801,p1,'r',1:1801,p2,'r')

 
load V_Rtot
load konRtot
load koff
load kdeg

xlabel('time')
ylabel('Order of accuracy, P')
title([' nu/ R_t_o_t = ',num2str(V_Rtot),'   k_o_nR_t_o_t = ',num2str(konRtot), '  k_o_f_f = ',num2str(koff)])
