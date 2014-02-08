function f = coor(aero)

m=aero(1);
p=aero(2);
t=aero(3);
n=aero(4);

%generating x coordinates using half cosine spacing
chordx=[];
for j=0:n
    x=(1-cos(j*pi/n))/2; %full cosine spacing
    %x=1-cos(j*pi/(2*n)); %half cosine spacing
    chordx=[chordx;x;];
end

%NACA Specific Code
%Code for Symmetric Airfoils
if(m==0&&p==0)
    y=(t/20)*((0.2969*(chordx.^0.5))-(.126*chordx)-(0.3516*(chordx.^2))+(.2843*(chordx.^3))-(0.1036*(chordx.^4)));
    y(n+1)=0;
    xupper=flipud(chordx);
    xlower=chordx;
    yupper=flipud(y);
    ylower=-y;
    xy=[xupper yupper;xlower ylower];
    xy(n+1,:)=[];
    f=xy;
    %plot(xy(:,1),xy(:,2))
    %axis equal
else
%Code for Cambered Airfoils    
%Camber Line
ycamb=[];
diffc=[];

for j=1:(n+1)
    if chordx(j) <= (p/10)
        x=m*chordx(j)*((p/5)-chordx(j))/(p^2);
        y=(m*((p/5)-chordx(j))/p^2)-(m*chordx(j)/(p^2));
        ycamb=[ycamb; x;];
        diffc=[diffc; y;];
    else
        x=m*(1-chordx(j))*(1+chordx(j)-p/5)/((10-p)^2);
        y=(m*(1-chordx(j))/((10-p)^2))-(m*(1-(p/5)+chordx(j))/((10-p)^2));
        ycamb=[ycamb; x;];
        diffc=[diffc; y;];
    end
end

theta=atan(diffc);
%Thickness about the mean line
ythick=(t/20)*((0.2969*(chordx.^0.5))-(.126*chordx)-(0.3516*(chordx.^2))+(.2843*(chordx.^3))-(0.1036*(chordx.^4)));
ythick(n+1)=0;

%GENERATE Airfoil Surface Coordinates
%upper surface
xupper=chordx-(ythick.*sin(theta));
yupper=ycamb+(ythick.*cos(theta));
%lower surface
xlower=chordx+(ythick.*sin(theta));
xlower=flipud(xlower);
ylower=ycamb-(ythick.*cos(theta));
ylower=flipud(ylower);
xy=[xlower ylower;xupper yupper];
xy(n+1,:)=[];
f = xy;
%plot(xy(:,1),xy(:,2))
%axis equal
end
%END of NACA Specific Module


