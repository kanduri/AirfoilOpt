%MAIN
clc
%airfoil parameters
aero=[2 4 12];

%no. of panels
number=10;
n=(number/2);

aero=[aero n];

xy=coor(aero);

%specify inlet conditions;
v=20; %velocity
alpha=4; %angle of attack
vel = [v*cos(alpha*pi/100) v*sin(alpha*pi/100)];

mid=[xy(:,1) xy(:,2)];
normals=[xy(:,3) xy(:,4)];
xcor=xy(:,5);
ycor=xy(:,6);

%declare a blank infuence matrix
inf=zeros(((2*n)+1),((2*n)+1));

%constant strength source/vortex distribution calculations

for j=1:2*n %panel at which velocity is being influenced
    for k=1:2*n %panel due to which velocity is being influenced at jth panel
        if j==k
            inf(j,k)=0.5;
        else
            theta=atan((ycor(k+1)-ycor(k))/(xcor(k+1)-xcor(k)));
            theta1=atan(mid(j,2)/(mid(j,1)-xcor(k)))-theta;
            theta2=atan(mid(j,2)/(mid(j,1)-xcor(k+1)))-theta;
            r1sq=(mid(j,1)-xcor(k))^2+(mid(j,2)-ycor(k))^2;
            r2sq=(mid(j,1)-xcor(k+1))^2+(mid(j,2)-ycor(k+1))^2;
            u=log(r1sq/r2sq)/(2*pi);
            w=(theta2-theta1)/(2*pi);
            inf(j,k)=[u w]*normals(j,:)';
            inf(j,(2*n)+1)=inf(j,n+1)+[w u]*normals(j,:)';
        end
    end
end
%last row calculations for kutta condition
for j=1:2*n
    inf((2*n)+1,j)=inf(j,1)+inf(j,2*n);
end
inf((2*n)+1,(2*n)+1)=inf(1,(2*n)+1)+inf((2*n),(2*n)+1);

%velocity component vector
V=zeros(2*n+1,1);
for j=1:2*n
    V(j)=-vel*normals(j,:)';
end

%source and vortex strengths
sigma=(inf^-1)*V;
sigma
%inf
%All source and vortex based influence coefficients are calculated.


%Function 

%COOR.M

function f = coor(aero)

m=aero(1);
p=aero(2);
t=aero(3);
n=aero(4);

%generating x coordinates
chordx=[];
for j=0:n
    x=(1-cos(j*pi/n))/2;
    chordx=[chordx;x;];
end
x=[];
y=[];
%NACA Specific Code
%Camber Line
ycamb=[];
diffc=[];
theta=[];
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
x=[];
y=[];
theta=atan(diffc);
%Thickness about the mean line
ythick=(t/20)*((0.2969*(chordx.^0.5))-(.126*chordx)-(0.3516*(chordx.^2))+(.2843*(chordx.^3))-(0.1036*(chordx.^4)));
ythick(n+1)=0;

%GENERATE Airfoil Surface Coordinates
%upper surface
xupper=chordx-(ythick.*sin(theta));
xupper=flipud(xupper);
yupper=ycamb+(ythick.*cos(theta));
yupper=flipud(yupper);
%lower surface
xlower=chordx+(ythick.*sin(theta));
ylower=ycamb-(ythick.*cos(theta));
xy=[xupper yupper;xlower ylower];
xy(n+1,:)=[];
%END of NACA Specific Module

%Following Code is general. Not specific to NACA module and must be shifted
%to influence module which is parametrization independent

%Generate Midpoints and Normal Unit Vectors
mid=[];
z=[];
slope=[];
for j=1:n
    x=(xy(j,1)+xy(j+1,1))/2;
    y=(xy(j,2)+xy(j+1,2))/2;
    z=(xy(j,2)-xy(j+1,2))/(xy(j,1)-xy(j+1,1));
    mid=[mid; x y];
    slope=[-slope;z];
end
for j=(n+1):(2*n)
    x=(xy(j,1)+xy(j+1,1))/2;
    y=(xy(j,2)+xy(j+1,2))/2;
    z=(xy(j,2)-xy(j+1,2))/(xy(j,1)-xy(j+1,1));
    mid=[mid; x y];
    slope=[slope;z];
end
nx=(slope./(((slope.^2)+1)).^0.5);
normal=[nx -nx./slope;0 0];
slope=[slope;0];
mid=[mid;0 0];
f=[mid normal xy slope];
