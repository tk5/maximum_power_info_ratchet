function [tau] = numerical_MFPT(V,xGrid,a,b)
%NUMERICAL_MFPT calculates the mean first passage time in the potential V
%specified on the grid xGrid from a to b numerically
%
% INPUTS: 
%       V: vector containing potential
%   xGrid: vector containing x-grid
%       a: origin
%       b: boundary
%
% OUTPUTS:  
%   tau: mean first passage time
%
% author:   JEhrich
% version:  1.3 (2021-03-29)
% changes:  fixed indexing in potential by swapping i_a and i_b

% calculate index of threshold
[~,i_a] = min(abs(a-xGrid));
% calculate index of reset
[~,i_b] = min(abs(b-xGrid));

% compute MFPT
dx = diff(xGrid(1:2));
% inner integral
ker = cumsum(exp(-V)*dx);
% outer integral
tau = sum(dx*exp(V(i_a:i_b)).*ker(i_a:i_b));



end