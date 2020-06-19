function [x,w] = GaussLegendre(N,a,b)
%% GaussLegendre Computes Legendre-Gauss nodes (x) and weights (w) on interval [a,b] with truncation order N.
% Suppose you have a continuous function f(x) defined on [a,b], that you can evaluate at any x in [a,b].
% Compute the definite integral using w'*f(x);
%
% INPUTS
% N    (1,1)  number of quadrature points
% a    (M,1)  lower bound of quadrature
% b    (M,1)  upepr bound of quadrature
% NOTE: either 'a' or 'b' must be a scalar
% 
% OUTPUTS
% x    (N,M)  vector of quadrature locations (in ascending order)
% w    (N,M)  vector of quadrature weights
% 
%
% Adapted from Greg von Winckel - 02/25/2004
% 
% AUTHORS: M. Iacopini, F. Ravazzolo, and L. Rossini 
% 
% TITLE: "Measuring and evaluating asymmetry in density forecasting"
% 
% AVAILABLE ON ....
% 
% PLEASE CITE AS: Iacopini,M., Ravazzolo, F. & Rossini, L. (2020) - "Measuring and evaluating asymmetry in density forecasting",
% available at .....
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N  = N-1;
N1 = N+1;
N2 = N+2;

xu = linspace(-1,1,N1)';
% Initial guess
y = cos((2*(0:N)'+1)*pi/(2*N+2)) + (0.27/N1)*sin(pi*xu*N/N2);
% Legendre-Gauss Vandermonde Matrix
L = zeros(N1,N2);
% Derivative of LGVM
Lp = zeros(N1,N2);
% Compute the zeros of the N+1 Legendre Polynomial using the recursion relation and the Newton-Raphson method
y0 = 2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
   L(:,1)  = 1;
   Lp(:,1) = 0;
   L(:,2)  = y;
   Lp(:,2) = 1;
   for k=2:N1
      L(:,k+1) = ( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) ) /k;
   end
   Lp = (N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
   y0 = y;
   y  = y0-L(:,N2)./Lp;
end
% Linear map from[-1,1] to [a,b]
x = (a.*(1-y)+b.*(1+y))/2;
% Compute the weights
w = (b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

% sort abscissas in ascending order (and adjust weigths accordingly) 
[x,idx] = sort(x,'ascend');
w = w(idx);
end
