function [x,y] = naca4(naca,npanel,x,y)


m=naca(1)/100;
p=naca(2)/10;
t=naca(3)/100;

% compute thickness and camber distributions

if mod(npanel,2) ~= 0
    sprintf('Please choose an even number of panels');
    sprintf('Exiting...');
    exit;
end

nside = npanel / 2 +1;

% camber distribution

for i=1:nside
    xx(i) = (1-cos(i*pi/nside))/2; %full cosine spacing
    %xx(i) = 1-cos(i*pi/(2*nside)); %half cosine spacing
    yt(i) = ( 0.29690*sqrt(xx(i)) -0.12600*xx(i)     ...
             -0.35160*xx(i)^2      +0.28430*xx(i)^3  ... 
             -0.10360*xx(i)^4) * t / 0.20;
    if xx(i) < p
        yc(i) = m/p^2 * (2*p*xx(i) -xx(i)^2);
    else
        yc(i) = m/(1 -p)^2 * ((1 -2*p) + 2*p*xx(i)-xx(i)^2);
    end
end

% airfoil shape = camber + thickness

for i=1:nside
    x(nside+i-1) = xx(i);
    x(nside-i+1) = xx(i);
    y(nside+i-1) = yc(i) +yt(i);
    y(nside-i+1) = yc(i) -yt(i);

end

return
