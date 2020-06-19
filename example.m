%% Code for examples of Asymmetric Continuous Probabilisty Score (ACPS)
%
% Written by 
% 
% AUTHORS: M. Iacopini, F. Ravazzolo, and L. Rossini 
% 
% TITLE: "Measuring and evaluating asymmetry in density forecasting"
% 
% AVAILABLE ON ....
% 
% PLEASE CITE AS: Iacopini,M., Ravazzolo, F. & Rossini, L. (2020) - "Measuring and evaluating asymmetry in density forecasting",
% available at .....
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
path = matlab.desktop.editor.getActiveFilename;
sf = strfind(path,'/');   addpath(genpath(path(1:sf(end))));   cd(path(1:sf(end)));
set(groot,'defaultAxesTickLabelInterpreter','latex','defaulttextinterpreter','latex','defaultLegendInterpreter','latex');
set(0,'DefaultFigurePosition',[100,500,400,400],'DefaultAxesFontSize',20,'DefaultTextFontSize',20);

dire = [cd,'/'];
res  = '-r120';
% information about random number generator -- to replicate the results
rng(20384);


%% (1) Target: NORMAL  --  Forec: NORMAL

% forecaster:  N(mu,s^2)
mu = [0.0, -3.0, 3.0, 0.0];
s  = [1.0,  1.0, 1.0, 4.0];

% generate observations
N   = 100;   % number of observations
M   = 100;   % number of values from forecasting distribution (to compute Empirical CDF)
mu0 = mu(1);     s0 = s(1);
obs = sort( mu0 + s0*randn([N,1]) );
% level of 'asymmetry':   c = 0.50 --> symmetric loss
c = linspace(0.05,0.95,5);

% compute score for each observation (the higher the better)
score = zeros(length(obs),length(c),length(mu));
for i=1:length(c)
   for j=1:length(mu)
      forec = mu(j) + s(j)*randn([M,1]);
      score(:,i,j) = ACPS(forec, obs, c(i));
   end
   disp(['Score ',num2str(i),' of ',num2str(length(c)),' computed']);
end
% compute SCORE (average over observations)
SA = squeeze(mean(score,1)) * 1e-3;
posMax = zeros(length(c),1);
SMax   = zeros(length(c),1);
for i=1:length(c)
   [SMax(i),posMax(i)] = max(SA(i,:));
end
% figure('position',[100,500,600,250]);
% plot(c,posMax,'-o','LineWidth',1.5); xlim([c(1),c(end)]); ylim([0,length(mu)+1]);

% show best model
disp(['Best model -->  ', num2str(posMax')]);
% collect all scores
[~,pos_best_NN] = sort(SA,2,'descend');


%%%% PLOTS %%%%
xx = linspace(-5,5,2000);   yy = linspace(-25,25,2000);
figure('Position',[100,500,500,300]); box on; hold on;
for j=1:length(mu)
   plot(yy,normcdf(yy,mu(j),s(j)),'LineWidth',1.6);
end
plot(xx,normcdf(xx,mu0,s0),'-k','LineWidth',1.6); xlim([-5,5]); ylim([0,1]);
print(gcf,[dire,'simulation_ACPS_NormNorm'],'-depsc2',res,'-opengl'); close(gcf);



%% (2) Target: STUDENT-t  --  Forec: STUDENT-t

% forecaster:  t(nu)
mu = [-3.0, 2.0, 0.0,  4.0];
s  = [ 1.0, 1.0, 1.0,  1.0];
nu = [ 3.0, 3.0, 5.0, 15.0];

% generate observations
N   = 100;   % number of observations/forecasts
M   = 200;   % number of values from forecasting distribution (to compute Empirical CDF)
mu0 = mu(3);   s0 = s(3);    nu0 = nu(3);
obs = sort( random('tLocationScale',mu0,s0,nu0,[N,1]) );
% level of 'asymmetry':   c = 0.50 --> symmetric loss
c = linspace(0.05,0.95,5);

score = zeros(length(obs),length(c),length(mu));
for i=1:length(c)
   for j=1:length(mu)
      forec = random('tLocationScale', mu(j), s(j), nu(j), [M,1]);
      score(:,i,j) = ACPS(forec, obs, c(i));
%       score(:,i,j) = S_t(obs, c(i), mu(j), s(j), nu(j));
   end
   disp(['Score ',num2str(i),' of ',num2str(length(c)),' computed']);
end
% compute SCORE (average over observations)
SA = squeeze(mean(score,1)) * 1e-3;
posMax = zeros(length(c),1);
SMax   = zeros(length(c),1);
for i=1:length(c)
   [SMax(i),posMax(i)] = max(SA(i,:));
end
% figure('position',[100,500,600,250]);
% plot(c,posMax,'-o','LineWidth',1.5); xlim([c(1),c(end)]); ylim([0,length(a)+1]);

% show best model
disp(['Best model -->  ', num2str(posMax')]);
% collect all scores
[~,pos_best_tt] = sort(SA,2,'descend');


%%%% PLOTS %%%%
xx = linspace(-20,20,2000);   yy = linspace(-20,20,2000);
figure('Position',[100,500,500,300]); box on; hold on;
for j=1:length(mu)
   plot(yy,cdf('tLocationScale',yy,mu(j),s(j),nu(j)),'LineWidth',1.6);
end
plot(xx,cdf('tLocationScale',xx,mu0,s0,nu0),'-k','LineWidth',1.6); xlim([-10,10]); ylim([0,1]);
print(gcf,[dire,'simulation_ACPS_tt'],'-depsc2',res,'-opengl'); close(gcf);




%% (3) Target: GAMMA  --  Forec: GAMMA

% forecaster:  Ga(a,b)
a = [1.0, 2.0, 1.5, 1.0];
b = [1.0, 1.0, 1.5, 2.0];

% generate observations
N   = 100;   % number of observations/forecasts
M   = 200;   % number of values from forecasting distribution (to compute Empirical CDF)
a0 = a(2);     b0 = s(2);
obs = sort( gamrnd(a0,b0,[N,1]) );
% level of 'asymmetry':   c = 0.50 --> symmetric loss
c = linspace(0.05,0.95,5);

score = zeros(length(obs),length(c),length(a));
for i=1:length(c)
   for j=1:length(a)
      forec = gamrnd(a(j),b(j),[M,1]);
      score(:,i,j) = ACPS(forec, obs, c(i));
   end
   disp(['Score ',num2str(i),' of ',num2str(length(c)),' computed']);
end
% compute SCORE (average over observations)
SA = squeeze(mean(score,1)) * 1e-3;
posMax = zeros(length(c),1);
SMax   = zeros(length(c),1);
for i=1:length(c)
   [SMax(i),posMax(i)] = max(SA(i,:));
end
% figure('position',[100,500,600,250]);
% plot(c,posMax,'-o','LineWidth',1.5); xlim([c(1),c(end)]); ylim([0,length(a)+1]);

% show best model
disp(['Best model -->  ', num2str(posMax')]);
% collect all scores
[~,pos_best_GG] = sort(SA,2,'descend');


%%%% PLOTS %%%%
xx = linspace(0,5,2000);   yy = linspace(0,25,2000);
figure('Position',[100,500,500,300]); box on; hold on;
for j=1:length(mu)
   plot(yy,gamcdf(yy,a(j),b(j)),'LineWidth',1.6);
end
plot(xx,gamcdf(xx,a0,b0),'-k','LineWidth',1.6); xlim([0,5]); ylim([0,1]);
print(gcf,[dire,'simulation_ACPS_GaGa'],'-depsc2',res,'-opengl'); close(gcf);



%% (4) Target: BETA  --  Forec: BETA

% forecaster:  Be(a,b)
a = [1.0, 1.0, 1.0, 5.0];
b = [1.0, 5.0, 2.0, 5.0];

% generate observations
N   = 100;   % number of observations
M   = 100;   % number of values from forecasting distribution (to compute Empirical CDF)
a0 = a(3);     b0 = b(3);
obs = sort( betarnd(a0,b0,[N,1]) );
% level of 'asymmetry':   c = 0.50 --> symmetric loss
c = linspace(0.05,0.95,5);

% compute score for each observation (the higher the better)
score = zeros(length(obs),length(c),length(a));
for i=1:length(c)
   for j=1:length(a)
      forec = betarnd(a(j),b(j),[M,1]);
      score(:,i,j) = ACPS(forec, obs, c(i));
   end
   disp(['Score ',num2str(i),' of ',num2str(length(c)),' computed']);
end
% compute SCORE (average over observations)
SA = squeeze(mean(score,1)) * 1e-3;
posMax = zeros(length(c),1);
SMax   = zeros(length(c),1);
for i=1:length(c)
   [SMax(i),posMax(i)] = max(SA(i,:));
end
% figure('position',[100,500,600,250]);
% plot(c,posMax,'-o','LineWidth',1.5); xlim([c(1),c(end)]); ylim([0,length(mu)+1]);

% show best model
disp(['Best model -->  ', num2str(posMax')]);
% collect all scores
[~,pos_best_BB] = sort(SA,2,'descend');


%%%% PLOTS %%%%
xx = linspace(0,1,2000);   yy = linspace(0,1,2000);
figure('Position',[100,500,500,300]); box on; hold on;
for j=1:length(a)
   plot(yy,betacdf(yy,a(j),b(j)),'LineWidth',1.6);
end
plot(xx,betacdf(xx,a0,b0),'-k','LineWidth',1.6); xlim([0,1]); ylim([0,1]);
print(gcf,[dire,'simulation_ACPS_BeBe'],'-depsc2',res,'-opengl'); close(gcf);

