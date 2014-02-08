clc
%airfoil parameters
aero=[4 4 12];

%no. of panels
number=10;
n=(number/2);

aero=[aero n];

xy=coor(aero);
xy=midnorm(xy,n);
%specify inlet conditions;
v=1; %velocity
alpha=10; %angle of attack
al=alpha*pi/180; %convert to radian

mid=[xy(:,1) xy(:,2)];
normals=[xy(:,3) xy(:,4)];
xcor=xy(:,5);
A=xcor;A(2*n+1)=[];
ycor=xy(:,6);

%define tangent vectors for each panel along the freestream velocity
tangents=fliplr(normals);
tangents=[-ones(n,1) ones(n,1); ones(n,1) -ones(n,1);0 0].*tangents;

%declare a blank infuence matrix
inf=zeros(((2*n)+1),((2*n)+1)); %normal velocity
inft=zeros((2*n),((2*n)+1));    %tangential velocity

%constant strength source/vortex distribution calculations
dummy=0;
st=zeros(2*n+1);
ct=zeros(2*n+1);

st=-normals(:,1);
ct=normals(:,2);

pi2inv = (1/(2*pi));
for j=1:2*n %panel at which velocity is being influenced
    for k=1:2*n %panel due to which velocity is being influenced at jth panel
        if j==k
            inf(j,k)=0.5;
            inft(j,2*n+1)=inft(j,2*n+1)+0.5;
        else
            
            ctdkj        = ct(j)*ct(k) +st(j)*st(k);
            stdkj        = st(j)*ct(k) -st(k)*ct(j);
            
            dxj    = mid(j,1)-xcor(k);
            dxjp   = mid(j,1)-xcor(k+1);
            dyj    = mid(j,2)-ycor(k);
            dyjp   = mid(j,2)-ycor(k+1);
            
 
            u=0.5*log((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj))*pi2inv; %log term
            w=atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj)*pi2inv; %tan term
            
            dummy1=[u w]*[stdkj ctdkj]';
            inf(j,k)=dummy1;
            dummy2=[u w]*[ctdkj -stdkj]';
            
            %last row calculations for kutta condition
            if j==1||j==2*n
                inf(2*n+1,k)=inf(2*n+1,k)-dummy2;
            end
            %last column calculations for vortex
            inf(j,2*n+1)=inf(j,2*n+1)+dummy2;
            %tangential velocity influence calculations
            inft(j,k)=-dummy2;
            inft(j,2*n+1)=inft(j,2*n+1)+dummy1;
        end
    end
end

%last term calculations for kutta condition
for j=1:2*n
    inf((2*n)+1,(2*n)+1)=inf((2*n)+1,(2*n)+1)+inf(1,j)+inf(2*n,j);
end

%velocity component vector
V=zeros(2*n+1,1);

for j=1:2*n
    V(j)=st(j)*cos(al) -sin(al)*ct(j);
end

V(2*n+1)=-(ct(1)     *cos(al) +st(1)     *sin(al)) ...
              -(ct(2*n)*cos(al) +st(2*n)*sin(al));


%source and vortex strengths
sigma=(inf^-1)*V;

%All source and vortex based influence coefficients are calculated.

%Calculation of Tangential Velocities at center of each panel
tangents(2*n+1,:)=[];
vt=(inft*sigma);
cp=1.-((vt./v).^2);
plot(A,-cp)