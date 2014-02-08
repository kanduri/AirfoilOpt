function vt = veldis(qg,x,y,xbar,ybar,st,ct,al,npanel)

% flow tangency boundary condition - source distribution

for i=1:npanel
    vt(i) = ct(i)*cos(al) +st(i)*sin(al);
    for j=1:npanel
        rp     = sqrt((xbar(i) -x(j+1))^2 +(ybar(i) -y(j+1))^2);
        rm     = sqrt((xbar(i) -x(j)  )^2 +(ybar(i) -y(j)  )^2);
        betaij = atan2 ((ybar(i) -y(j+1))*(xbar(i) -x(j)) ...
                -(xbar(i) -x(j+1))*(ybar(i) -y(j)),  ...
                    (xbar(i) -x(j+1))*(xbar(i) -x(j))   ... 
                       +(ybar(i) -y(j+1))*(ybar(i) -y(j)));
        
                   %calculation of tangent influence coefficients
        if i==j
            log_term  = 0.0;
            beta_term = pi;
        else
            log_term  = log(rp/rm);
            beta_term = betaij;
        end
        
        ctimtj = ct(i)*ct(j) +st(i)*st(j);
        stimtj = st(i)*ct(j) -st(j)*ct(i);
        
        vt(i) = vt(i) + (qg(j)/(2*pi))*(stimtj*beta_term -ctimtj*log_term) ...
                + (qg(npanel+1)/(2*pi))*(stimtj*log_term +ctimtj*beta_term);
    end
end

return