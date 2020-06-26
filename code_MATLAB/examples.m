%% test code for plotting an Half-Disk Denisty Strip (HDDS)
clear;
path = matlab.desktop.editor.getActiveFilename;  sf=strfind(path,'/');
addpath(genpath(path(1:sf(end)))); cd(path(1:sf(end)));
set(groot,'defaultAxesTickLabelInterpreter','latex','defaulttextinterpreter','latex','defaultLegendInterpreter','latex');
set(0,'DefaultFigurePosition',[100,500,400,400],'DefaultAxesFontSize',20,'DefaultTextFontSize',20);


%% example 1 - Normal
% setup
diameter = 1.00;        % diameter of HDDS
radius = diameter / 2;
Nbins  = 100;           % number of bins of the HDDS
Ndata  = 500;           % number of simulated data points
inv_scale = 4.0;        % scaling factor for color shading

% generate data
x = 2.0*randn([Ndata,1]);
discrete = false;
% specify bounds of HDDS -- use NaN(2,1) for default  bounds = [min(x),max(x)];
bounds = [-4,4];
% select color palette
ramp = flipud(gray(300));
legend = true;

% plot HDDS (get function handle)
h1 = HDDS(x,bounds,discrete,Nbins,ramp,inv_scale,radius,legend);



%% example 2 - Gamma
% setup
diameter = 1.00;        % diameter of HDDS
radius = diameter / 2;
Nbins  = 200;           % number of bins of the HDDS
Ndata  = 500;           % number of simulated data points
inv_scale = 4.0;        % scaling factor for color shading

% generate data
x = gamrnd(3.0, 2.0, [Ndata,1]);
discrete = false;
% specify bounds of HDDS -- use NaN(2,1) for default  bounds = [min(x),max(x)];
bounds = NaN;
% select color palette
ramp = flipud(hot(300));
legend = true;

% plot HDDS (get function handle)
h2 = HDDS(x,bounds,discrete,Nbins,ramp,inv_scale,radius,legend);



%% example 3 - Negative Binomial
% setup
diameter = 1.00;        % diameter of HDDS
radius = diameter / 2;
Nbins  = 100;           % number of bins of the HDDS
Ndata  = 500;           % number of simulated data points
inv_scale = 4.0;        % scaling factor for color shading

% generate data
x = nbinrnd(3.0, 0.4, [Ndata,1]);
discrete = true;
% specify bounds of HDDS -- use NaN(2,1) for default  bounds = [min(x),max(x)];
bounds = NaN;
% select color palette
ramp = flipud(gray(300));
legend = true;

% plot HDDS (get function handle)
h3 = HDDS(x,bounds,discrete,Nbins,ramp,inv_scale,radius,legend);

