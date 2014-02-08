function ainfl = infl_coeff(x,y,xbar,ybar,st,ct,ainfl,npanel)

pi2inv = 1 / (2*pi);

% set vn = 0 at midpoint of ith panel

for i=1:npanel

    % find contribution of the jth panel
    
    for j=1:npanel
        
        log_term = 0.0;
        tan_term = pi;

        if i ~= j
            dxj    = xbar(i) -x(j);
            dxjp   = xbar(i) -x(j+1);
            dyj    = ybar(i) -y(j);
            dyjp   = ybar(i) -y(j+1);
            
            log_term = 0.5*log((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj));
            tan_term = atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj);
            
        end
        
        ctimtj        = ct(i)*ct(j) +st(i)*st(j);
        stimtj        = st(i)*ct(j) -st(j)*ct(i);
        ainfl(i,j)    = pi2inv*(tan_term*ctimtj +log_term*stimtj);
        add           = pi2inv*(log_term*ctimtj -tan_term*stimtj);
        ainfl(i,npanel+1) = ainfl(i,npanel+1) + add; %last column
        
        if (i == 1 | i == npanel)

            % if the ith panel touches trailing edge
            %  add contribution to kutta condition
            
            ainfl(npanel+1,j)        = ainfl(npanel+1,j) - add;
            ainfl(npanel+1,npanel+1) = ainfl(npanel+1,npanel+1) + ainfl(i,j);
        end
    end
end 
 
if rank(ainfl) ~= npanel+1
    exit
end

return
