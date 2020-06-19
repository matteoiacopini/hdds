function [ACPS,ACPS_mean] = ACPS(forec, obs, c)
%% ACPS Computes ACPS for each observation
% INPUTS
%  forec  (M,1)  draws from forecasting distribution (which has to be estimated via empirical CDF)
%  obs    (N,1)  observed values
%  c      (1,1)  asymmetry level in (0,1), where:
%                  c = 0.5 means symmetric loss,
%                  c < 0.5 penalises right-shifted CDF,
%                  c > 0.5 penalises left-shifted CDF
% 
% OUTPUTS
%  ACPS       (N,1) ACPS at each observation
%  ACPS_mean  (1,1) ACPS (mean across all observations)
%
% Written by
% AUTHORS: M. Iacopini, F. Ravazzolo, and L. Rossini 
% 
% TITLE: "Measuring and evaluating asymmetry in density forecasting"
% 
% AVAILABLE ON ....
% 
% PLEASE CITE AS: Iacopini,M., Ravazzolo, F. & Rossini, L. (2020) - "Measuring and evaluating asymmetry in density forecasting",
% available at .....
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isvector(obs) && size(obs,1)~=1
   obs = obs';
end
N = length(obs);
ACPS = zeros(N,1);

% info for quadrature
Nquad = 500;      % number of points
yu    =  1000;    % upper bound (truncation)
yl    = -1000;    % lower bound (truncation)

% generate Gauss quadrature locations and weights
[xx_lo,ww_lo] = GaussLegendre(Nquad,yl,obs);    % Nquad x N --> quadrature for each value of 'obs'
[xx_up,ww_up] = GaussLegendre(Nquad,obs,yu);

% compute empirical CDF of forecasts
[ff,zz] = ecdf(forec);
FF = repmat(ff,[1,Nquad]);

% compute Score for each observation
for ii=1:N   
   % compute Empirical CDF at 'xx'
   Fx_lo = max(FF.*(zz <= xx_lo(:,ii)'),[],1);
   Fx_up = max(FF.*(zz <= xx_up(:,ii)'),[],1);
   
   % compute 'right' integral on (obs,+inf)
   tmp_up = ((1-c)^2 - (1-Fx_up).^2) .* ((1-c)^(-2)*(Fx_up > c) + c^(-2)*(Fx_up <= c));
   int_up = sum(ww_up(:,ii).*tmp_up', 1);
   % compute 'left'  integral on (-inf,obs)
   tmp_lo = (c^2 - Fx_lo.^2) .* ((1-c)^(-2)*(Fx_lo > c) + c^(-2)*(Fx_lo <= c));
   int_lo = sum(ww_lo(:,ii).*tmp_lo', 1);

   ACPS(ii,1) = int_lo + int_up;
end

ACPS_mean = mean(ACPS);
end