% Files:
%     - experiments/test_Ps.m: Calculates the P matrix and GP gain using 
%         different methods. We can see that the results from different
%         methods are consistant however with a small r (e.g. r=10), GP does 
%         not produce GLQR. But as r increases (e.g. r=100) GP approaches
%         GLQR.
%     - experiments/test_Ss.m: Same as experiments/test_Ps.m with similar
%         results but for S matrix and GL gain. Notice that in this case we
%         cannot increase the r too much because the number of parameters
%         (size of the S matrix) will increase dramatically as r increases. 
%     - experiments/Eq58_Eq61_experiment.m: Shows that by changing RHS from
%         "RHS(i,:) = U_k" to "RHS(i,:) = U_k - (gamma^r)*U_kpr" we can
%         recreate the same exact result as analytical P from Numerical P.
%         This is an important result because it shows that Eq 30 is exact
%         and the extra part "(gamma^r)*U_kpr" only vanishes for the case
%         with large r and gamma<1. making gamma<1 a necessity if we are to
%         ignore the extra term (as they do in all RL literature). 
%         Yet, adding the extra term still does not make the results of
%         numerical P consistent for complex problems :(
%
%
%
%
%
%
