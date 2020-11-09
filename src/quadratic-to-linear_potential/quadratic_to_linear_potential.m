%QUADRATIC_TO_LINEAR_POTENTIAL produces evaluates output power for
%qaudratic-to-linar potential numerically. It then produces a figure
%illustrating the shape of the potential and plots the results for the
%output power
%
% OUTPUTS:  
%       1) creates figure with the output plots
%       2) creates txt file with plot data
%
% author:   JEhrich
% version:  1.2 (2020-11-04)
% changes:  added plot data output
clear 'all'
close 'all'
clc
% set interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% numerical parameters
% set size of grid
N = 1E6;
xGrid = linspace(-15,3,N);

%% plotting parameters
lW = 1.7; % linewidth
fS = 17; % font size

%% system parameters
% specify which xT values to evaluate
xTVec = logspace(-1.4,log10(max(xGrid)),200)';

% specify which f_1, f_2 combinations to probe
f1Vec = [1 1 1 1.3 1.3 1.3 2.8 2.8 2.8 inf];
f2Vec = [1 0.1 0.01 1 0.1 0.01 1 0.1 0.01 inf];
dg = 0.8;


%% A) example plot of potential
% x- axis
x_a = linspace(-5,5,1000);
% specify maximum slopes
f1 = 2;
f2 = 1;
% stitch together trapping potential
V0 = x_a.^2/2;
Vm = -f1*x_a - f1^2/2;
Vp = f2*x_a - f2^2/2;
Vt = nan(size(x_a));
Vt(x_a < -f1) = Vm(x_a < -f1);
Vt(-f1 <= x_a & x_a < f2) = V0(-f1 <= x_a & x_a < f2);
Vt(f2 <= x_a) = Vp(f2 <= x_a);

% create figure
figure('Position',[1,1,1000,450]);

% create axes for subfigure A
axes('Position',[0.05 0.12 0.4 0.78]);

% plot quadratic potential
plot(x_a,V0,'k--','linewidth',lW);
hold on;
% plot complete trap potential
plot(x_a,Vt,'k','linewidth',lW);
% plot markers where potential changes
plot(-f1, f1^2/2, '.k', 'MarkerSize',33)
plot(f2, f2^2/2, '.k', 'MarkerSize',33)

% set axis limits
axis([-4,4,0,7]);
% set axes font size
set(gca,'FontSize',fS);
% label x axis
xlabel('$x$','Interpreter','latex');
% add plot labels
text(3.2,4,'$\frac{1}{2}x^2$','FontSize',fS,'Interpreter','latex','horizontalAlignment','center');
text(3.2,2.1,'$\sim f_2\,x$','FontSize',fS,'Interpreter','latex','horizontalAlignment','center');
text(-3.3,3.2,'$\sim -f_1\,x$','FontSize',fS,'Interpreter','latex','horizontalAlignment','center');
% add figure label
text(-4.5,7.5,'A','FontSize',fS+9,'Interpreter','latex','horizontalAlignment','center');


%% compute power
% create empty vectors for outputs
tauVec = nan(length(xTVec),length(f1Vec));
xRVec = nan(length(xTVec),length(f1Vec)); % analytical reset position

% main loop
tic
for jj = 1:length(f1Vec)
    % output progress
    round(jj/length(f1Vec)*100)
    % set f1 and f2
    f1 = f1Vec(jj);
    f2 = f2Vec(jj);
    
    % loop through all thresholds
    for ii = 1:length(xTVec)
        xT = xTVec(ii);
        
        % compute analytic reset position
        if f1 >= f2 % case left side stronger than right
            if xT < f2
                xRVec(ii,jj) = xT;
            elseif f2 <= xT && xT < (f1^2+f2^2)/2/f2
                xRVec(ii,jj) = sqrt(2*f2*xT-f2^2);
            else
                xRVec(ii,jj) = (f1^2-f2^2+2*f2*xT)/2/f1;
            end
        else % case right side stronger than left
            if xT < f1
                xRVec(ii,jj) = xT;
            elseif f1 <= xT && xT < f2^2/2
                xRVec(ii,jj) = (f1^2+xT^2)/2/f1;
            else
                xRVec(ii,jj) = (f1^2 - f2^2 + 2*f2*xT)/2/f1;
            end
        end
        
        % parts of the trapping potential
        V0 = xGrid.^2/2;
        Vm = -f1*xGrid - f1^2/2;
        Vp = f2*xGrid - f2^2/2;
        % stitch together trapping potential
        Vt = nan(size(xGrid));
        Vt(xGrid < -f1) = Vm(xGrid < -f1);
        Vt(-f1 <= xGrid & xGrid < f2) = V0(-f1 <= xGrid & xGrid < f2);
        Vt(f2 <= xGrid) = Vp(f2 <= xGrid);
        % complete potential
        V = Vt + dg*xGrid;
        
        % compute MFPT
        tauVec(ii,jj) = numerical_MFPT(V,xGrid,xT,-xRVec(ii,jj)); 
    end
end
toc

% compute power
PVec = dg*(xTVec+xRVec)./tauVec;
    


%% B) plot output power
% specify linestyles
linestyle{1} = '--';
linestyle{2} = '-.';
linestyle{3} = ':';
linestyle{4} = '--';
linestyle{5} = '-.';
linestyle{6} = ':';
linestyle{7} = '--';
linestyle{8} = '-.';
linestyle{9} = ':';
% specify colors
colors = 'gggbbbrrr';

% create axes for subplot
axes('Position',[0.55 0.12 0.4 0.78]);

% plot curve for pure quadratic potential
semilogx(xTVec,PVec(:,end),'-k','linewidth',2.7);
hold on;
% plot invisible curves for legend
plot(nan,nan,'k--','linewidth',lW)
plot(nan,nan,'k-.','linewidth',lW)
plot(nan,nan,'k:','linewidth',lW)

% plot power for different f_1, f_2 combinations
for ii = 1:length(f1Vec)-1
    semilogx(xTVec,PVec(:,ii),'linestyle',linestyle{ii},'color',colors(ii),'linewidth',lW);
end
 
% add plot labels
text(0.1,0.268,'$f_1 = 2.8$','FontSize',fS,'Interpreter','latex','horizontalAlignment','center');
text(0.1,0.178,'$f_1 = 1.3$','FontSize',fS,'Interpreter','latex','horizontalAlignment','center');
text(0.1,0.09,'$f_1 = 1$','FontSize',fS,'Interpreter','latex','horizontalAlignment','center');
% add legend
legend({'$V_\mathrm{t}(x) = \frac{1}{2}x^2$','$f_2=1$','$f_2=0.1$','$f_2=0.01$'},'location','southwest');
% set font size
set(gca,'FontSize',fS);
% add axes labels
xlabel('$X_\mathrm{T}$','interpreter','latex')
ylabel('$P(k_\mathrm{B} T/\tau_\mathrm{r})$','interpreter','latex')
% set axis limits
axis([min(xTVec),max(xTVec),0,0.3]);
% add figure label
text(0.025,0.322,'B','FontSize',fS+9,'Interpreter','latex','horizontalAlignment','center');

% save figure
saveas(gcf, '../../doc/quadratic_to_linear','epsc')

%% write out data
fileID = fopen(['../../data/quadratic-to-linear_dg_' num2str(dg) '.txt'],'w');
% write comment line
fprintf(fileID,'%15s ','XT');
for ii = 1:length(f1Vec)
    fprintf(fileID,'%15s ',['f1=' num2str(f1Vec(ii)) ',f2=' num2str(f2Vec(ii))]);
end
fprintf(fileID,'\n');
% write out data
for ii = 1:length(xTVec)
    % write XT
    fprintf(fileID, '%10.9E ', xTVec(ii));
    % write P
    for jj = 1:length(f1Vec)-1
        fprintf(fileID, '%10.9E ', PVec(ii,jj));
    end
    fprintf(fileID, '%10.9E', PVec(ii,end));
    if ii < length(xTVec)
        fprintf(fileID,'\n');
    end
end
fclose(fileID);

