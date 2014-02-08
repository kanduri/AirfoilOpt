clc
angle=zeros(7,1);
clift=zeros(7,1);
naca=[4 2 12];
for j=1:7
    angle(j)=(j-1)*2.5;
    alpha=(j-1)*2.5;
    clift(j)=coeflift(naca,alpha);
end

plot(angle,clift)
axis equal