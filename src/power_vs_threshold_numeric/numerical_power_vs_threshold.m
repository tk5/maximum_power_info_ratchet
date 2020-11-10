%NUMERICAL_POWER_VS_THRESHOLD computes output power as afunction of
%threshold position for given gravitational offset dg
%
% OUTPUTS:  
%       creates .txt with tab delimiters of xT, power, velocity and MFPT
%
% author:   JEhrich
% version:  1.3 (2020-11-09)
% changes:  changed output file location to /data

clear 'all'
close 'all'
clc

% specify delta_g
dg = 0.8;

% specify values for XT
XTVec = linspace(1E-2,4,400)';

% emtpy vectors for outputs
tauVec = nan(size(XTVec));
vVec = nan(size(XTVec));
PVec = nan(size(XTVec));

tic
% main loop
for ii = 1:length(XTVec)
    ii
    XT = XTVec(ii);
    % grid of x values
    xGrid = linspace(-10,XT,1E6);
    % grid spacing
    dx = diff(xGrid(1:2));
    % potential
    V = xGrid.^2/2+dg*xGrid;
    % evaluate inner integral
    I_in = cumsum(dx*exp(-V));
    % find index of reset position -XT
    indMXT = find(xGrid<-XT,1,'last');
    % compute outer integral
    tauVec(ii) = sum(dx*exp(V(indMXT:end)).*I_in(indMXT:end));
    % evaluate velocity and power
    vVec(ii) = 2*XT/tauVec(ii);
    PVec(ii) = dg*vVec(ii);
end
toc

%% plot output power
figure();
plot(XTVec,PVec);

%% write out output power
%writematrix([XTs, P, v, tMFP],['num_power_threshold_dg_' num2str(dg) '.txt'],'Delimiter','tab')
fileID = fopen(['../../data/num_power_threshold_dg_' num2str(dg) '.txt'],'w');
% write comment line
fprintf(fileID,'%15s %15s %15s %15s\n','XT','P','v','tau');
% write out data
for ii = 1:length(XTVec)
    fprintf(fileID, '%10.9E %10.9E %10.9E %10.9E', XTVec(ii), PVec(ii), vVec(ii), tauVec(ii));
    if ii < length(XTVec)
        fprintf(fileID,'\n');
    end
end
fclose(fileID);

type(['../../data/num_power_threshold_dg_' num2str(dg) '.txt'])
