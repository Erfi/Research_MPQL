% This experiment if for drawing a heatmap for different
% values of r and gamma and checking the stability of the system. 
% -------------------------------------------------------------------------
clear all;
close all;

% --- System ---
[A,B,Cc,Dc,Q,R,Ac,Bc] = getSystemModel(4);

load('trussmodelforErfan.mat','dt');
[n,m] = size(B);
% --------------

r_range = [1,20:20:80];
gamma_range = [0.0:0.1:1.0];
stabilityMap = ones(length(r_range), length(gamma_range)) * -1;

if (1)
    
r_indx = 1;
for r=r_range
    gamma_indx = 1;
    for gamma=gamma_range
        r
        gamma
        S = calculateAnalyticalS(A,B,r,gamma,Q,R);
        [~,GL,~] = extractGainFromS(S,n,m);
        stabilityMap(r_indx,gamma_indx) =  checkStability(GL, A, B);
        
        %--- update indecies ---
        gamma_indx = gamma_indx + 1;
    end
    r_indx = r_indx + 1;
end

end
%% plotting
mymap = [1 0.2 0.2
    0.2 1 0.2];
 colormap(mymap);   % set colormap
 imagesc(gamma_range, r_range, stabilityMap);
 title('Stability Map for initial gain as r and gamma values change')
 xlabel('gamma')
 ylabel('r')
 set(gca,'xaxisLocation','top')
 grid on
 colorbar;          % show color scale
 
 %% plotting differently
 mymap = [1 0.2 0.2
    0.2 1 0.2];
h = heatmap(gamma_range, r_range,stabilityMap);
h.Colormap = mymap;
 
 
 title('Stability Map for initial gain as r and gamma values change')
 xlabel('gamma')
 ylabel('r')
%  set(gca,'xaxisLocation','top')
 grid on
 colorbar;          % show color scale
%  