function [l,st,ct,xbar,ybar] = panel_geometry(x,y,npanel)

% compute various geometrical quantities
%[x y]

for i=1:npanel
    l   (i) = sqrt((x(i+1) -x(i))^2 +(y(i+1) -y(i))^2);
    st  (i) = (y(i+1) -y(i))/l(i);
    ct  (i) = (x(i+1) -x(i))/l(i);
    xbar(i) = (x(i+1) +x(i))/2;
    ybar(i) = (y(i+1) +y(i))/2;
end

return




