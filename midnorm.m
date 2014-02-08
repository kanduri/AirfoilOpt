function f = midnorm(xy,n)
%Following Code is general. Not specific to NACA module. Parametrization
%independent.

%Generate Midpoints and Normal Unit Vectors
mid=[];
z=[];
sine=[];
for j=1:(2*n)
    x=(xy(j,1)+xy(j+1,1))/2;
    y=(xy(j,2)+xy(j+1,2))/2;
    z=(xy(j,2)-xy(j+1,2))/sqrt((xy(j,2)-xy(j+1,2))^2+(xy(j,1)-xy(j+1,1))^2);
    mid=[mid; x y];
    sine=[sine;z];
end
cosine=[-ones(n,1); ones(n,1)].*(((1-(sine.^2))).^0.5);
normal=[sine cosine;0 0];
mid=[mid;0 0];
f=[mid normal xy];