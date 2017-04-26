function [ x_kp1, done ] = cartpoleStep( x_k, u_k )
% This function takes action u_k in state x_k 
% and returns the resulting next state.
%
% Args:
%   x_k: Vector of current state 
%   u_k: Scalar value of input for the system
%
% Returns:
%   x_kp1: Vector of next time step state
    
    %constants
    gravity = 9.8;
    masscart = 1.0;
    masspole = 0.1;
    total_mass = (masscart + masspole);
    length = 0.5; %actually half the pole's length
    polemass_length = (masspole*length);
    force_mag = 1; 
    tau = 0.02; %seconds between state updates
    pos_threshold = 10;
    theta_threshold_radians = 15 * 2 * pi / 360;
    

    pos = x_k(1);
    pos_dot = x_k(2);
    theta = x_k(3);
    theta_dot = x_k(4);
    
    
    force = force_mag * u_k;
%     if (u_k == 1)  
%         force = force_mag;
%     else
%         force = -force_mag;
%     end
    
    costheta = cos(theta);
    sintheta = sin(theta);
    
    temp = (force + polemass_length * theta_dot * theta_dot * sintheta) / total_mass;
    thetaacc = (gravity * sintheta - costheta * temp) / (length * (4.0/3.0 - masspole * costheta * costheta / total_mass));
    posacc = temp - polemass_length * thetaacc * costheta / total_mass;
    pos = pos + tau * pos_dot;
    pos_dot = pos_dot + tau * posacc;
    theta = theta + tau * theta_dot;
    theta_dot = theta_dot + tau * thetaacc;
    
    x_kp1 = [pos, pos_dot, theta, theta_dot]';
    
    if ((pos < -pos_threshold) || (pos > pos_threshold) || (theta < -theta_threshold_radians) || (theta > theta_threshold_radians))
        done = 1;
    else
        done = 0;
    end
end

