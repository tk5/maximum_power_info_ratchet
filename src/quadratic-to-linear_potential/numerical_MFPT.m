function [tau] = numerical_MFPT(V,xGrid,a,b)
%NUMERICAL_MFPT calculates the mean first passage time in the potential V
%specified on the grid xGrid from a to b numerically
%
% INPUTS: 
%   V: vector containing potential
%   x: vector containing x-grid
%   a: origin
%   b: boundary
%
% OUTPUTS:  
%   tau: mean first passage time
%
% author:   JEhrich
% version:  1.1 (2020-11-04)
% changes:  code comments added, removed numerical calculation of reset
% positon

% calculate index of threshold
[~,i_a] = min(abs(a-xGrid));
% calculate index of reset
[~,i_b] = min(abs(b-xGrid));

% compute MFPT
dx = diff(xGrid(1:2));
% inner integral
ker = cumsum(exp(-V)*dx);
% outer integral
tau = sum(dx*exp(V(i_b:i_a)).*ker(i_b:i_a));



end

