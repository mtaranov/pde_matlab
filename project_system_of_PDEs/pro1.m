%%%%%%%%%%%%%% P R O J E C T  %%%%%%%%%%%%%%%%%
%1.- WRITE THE TRIDIAGONAL SOLVER
%2.- TEST TRIDIAGONAL SOLVER
%3.- WRITE PDE METHOD
%4.- TEST PDE METHOD IN LINEAR CASE
%5.- TEST METHOD ON PATTERN CASE
%6.- TEST DIFFERENT INITIAL CONDITIONS.

close all;
clear all;

%%%%%%%%%%%%   COEFFICIENTS AND BOUNDARY CONDITIONS  %%%%%%%%%%%%%%
a=1;
b=3;
uc=1;

n=1; 
dt = 0.05; 
n_end = 1/dt; % # of dt steps to take (time steps)

dx = 0.5;
j_end = 10/dx; % # of dx steps to move (space steps)
N=j_end;


Du=1; % diffusion coefficients
Dv=10;

%sig=0.015;
%rho=0.01;
%rho0=0.001;

mu1=Du*dt/dx^2;
mu2=Dv*dt/dx^2;
M=zeros(N,N);

% GENERATING THE 2 TRIDIAGONAL MATRICES
m=(N-1)/ 2;
Mu = (mu1+1)*eye(N) + diag((-mu1/2)*ones(2*m,1),1)+diag((-mu1/2)*ones(2*m,1),-1);
Mu(1,2) = -mu1;
Mu(N,N-1) = -mu1;

Mv = (mu2+1)*eye(N) + diag((-mu2/2)*ones(2*m,1),1)+diag((-mu2/2)*ones(2*m,1),-1);
Mv(1,2) = -mu2;
Mv(N,N-1) = -mu2;


%%CONSTRUCTING right hand side and then solving using thomas algorithm%%
%% Solving Mu * u = du (where du depends on gu)
%% Solving Mv * v = dv (where du depends on gu)

%%Initial u(n)
initial_u=1;
initial_v=10;

u = zeros(j_end,n_end);
u(:,1) = initial_u;

v = zeros(j_end,n_end);
v(:,1) = initial_v;

func1 = inline('a - b*u + cu/v ') %func1(a, b, uc, u,v);
func2 = inline('uc - v') %func2(b,u,v)


for n=1:n_end-1

    %Du = D*mu1;
    gu(1,n) = (1 - mu1)*u(1,n) + mu1*u(2,n);
    gv(1,n) = (1 - mu2)*v(1,n) + mu2*v(2,n);
    
    j=1; %this is the first space step
    %du(j,n) = 0;
    %dv(j,n) = 0;
    du(j,n) = gu(j,n) + dt*func1(a,b,uc, u(j,n), v(j,n));
    dv(j,n) = gv(j,n) + dt*func2(uc, v(j,n));
    
    j=j_end;
    gu(j,n) = (1 - mu1)*u(j,n) + mu1*u(j-1,n);
    du(j,n) = gu(j,n) + dt*func1(a,b,uc,(3*u(j,n)-u(j-1,n))/2 , (3*v(j,n) - v(j-1,n))/2);
    
    gv(j,n) = (1 - mu2)*v(j,n) + mu2*v(j-1,n);
    dv(j,n) = gv(j,n) + dt*func2(uc, (3*v(j,n) - v(j-1,n))/2);
    
    
    for j=2:j_end-1
        gu(j,n) = (1-mu1)*u(j,n) + mu1*(u(j+1,n) + u(j-1,n))/2;
        du(j,n) = gu(j,n) + dt*func1(a,b,uc,(3*u(j,n)-u(j-1,n))/2 , (3*v(j,n) - v(j-1,n))/2);
        %du(j,n) = gu(j,n) + dt*func1(a,(3*u(j,n)-u(j-1,n))/2 , (3*v(j,n) - v(j-1,n))/2);
        
        gv(j,n) = (1-mu2)*v(j,n) + mu2*(v(j+1,n) + v(j-1,n))/2;
        dv(j,n) = gv(j,n) + dt*func2(uc, (3*v(j,n) - v(j-1,n))/2);
        %dv(j,n) = gv(j,n) + dt*func2(b, (3*u(j,n)-u(j-1,n))/2, (3*v(j,n) - v(j-1,n))/2);
    end
    
    
    
    j=j+1; %this is the last space step
    
    %%%%%%%%%%%%%%%%%%%% THOMAS ALGORITHM (TRIDIAGONAL
    %%%%%%%%%%%%%%%%%%%% SOLVER)%%%%%%%%%%%%%%%
    %%%Solve Mu*u(:,n+1) = d(:,n)%%%
    temp_u = thomas(Mu,du(:,n));

    %temp_u = inv(Mu)*du(:,n);
    u(:,n+1) = temp_u;
    temp_v = thomas(Mv,dv(:,n));
    %temp_v = inv(Mv)*dv(:,n);
    v(:,n+1) = temp_v;
    
end

%u
%v

figure(1);
subplot(2,1,1)

mesh(u);
subplot(2,1,2)

mesh(v);


