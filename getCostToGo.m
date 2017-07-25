function [CTG, U] = getCostToGo( x_hist, u_hist, P, Q, R)
% plots the cost to go function over time
% using U(k) = x(k)'Qx(k) + u(k)'Ru(k) as utility/cost function.
if (P == false)
    len = min(size(x_hist,2),size(u_hist,2)); 
    CTG = nan;
    U = zeros(1,len); %utility
    for i=1:len
        U(i) = x_hist(:,i)'*Q*x_hist(:,i) + u_hist(:,i)'*R*u_hist(:,i);
    end    
else
    len = min(size(x_hist,2),size(u_hist,2)); 
    xu = vertcat(x_hist(:,1:len), u_hist(:,1:len));
    CTG = zeros(1,len); %cost to go
    U = zeros(1,len); %utility
    for i=1:len
        CTG(i) = xu(:,i)'*P*xu(:,i);
        U(i) = x_hist(:,i)'*Q*x_hist(:,i) + u_hist(:,i)'*R*u_hist(:,i);
    end
end

