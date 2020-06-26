function [col] = prob_to_col(dens,ramp,inv_scale)
%% prob_to_col maps probabilities to colors of a specified palette
   Ncol = size(ramp,1);
   [~,idx] = max( dens <= repmat(linspace(0,1/inv_scale,Ncol+1),[length(dens),1]),[],2);
   % corrections for max color
   idx(dens >= 1/inv_scale) = Ncol;
   idx(idx  >= Ncol) = Ncol;
   col = ramp(idx,:);
end