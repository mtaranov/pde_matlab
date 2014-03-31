clean;
to=0;
n_curves=20;%number of successive curves
t_curves=600;%time btw succesive curves; A-300, B-600, C,D-1800
tmax=t_curves*n_curves;%final time; 

ne=2;
n=100;%number of space steps
save('n','n')
co=zeros(1,ne*n);

[T,Y] = ode23t(@y_prime,to:t_curves:tmax,co);

for i=2:length(T)
    plot(100/(n-1)*(0:n-1),Y(i,2:2:n*ne))
    hold on
end
hold off

xlabel('distance (microns)')
ylabel('bound / R_t_o_t')


