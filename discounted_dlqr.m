function [ gain ] = discounted_dlqr( A,B,Q,R,N,gamma)
%This function uses dynamic programming to calculate the gain for 
%discounted LQR by solving the Algebraic Recatti Equation at each time 
%step backwards.
%
% Args:
%   A: System Dynamic Matrix
%   B: System Input Matrix
%   Q: State Weight Matrix
%   R: Input Weight Matrix
%   N: Number of steps until convergance. This needs to be large to
%      simulate infinite-horizon cost function
%   gamma: Discount or Forget Factor
%
% Returns:
%   gain: Gain for discounted discrete LQR
%--------------------------------------------------------------------------
    P_tp1 = (gamma^(N+1))*Q;
    for t=N:-1:0
        K_t = -pinv((gamma^(t))*R + B'*P_tp1*B)*B'*P_tp1*A;
        P_t = (gamma^(t))*Q + A'*P_tp1*A - A'*P_tp1*B*pinv((gamma^(t))*R + B'*P_tp1*B)*B'*P_tp1*A;
        P_tp1 = P_t;
    end
    gain = K_t;
end

