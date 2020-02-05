%This code is deisgned to plot the bifurcation diagram for the
%determinsitic version of the locust switching system in lecture 6
%By Kit Yates
%created 12/10/15
%Last modified 12/10/15

close all
clear all

%Define the rate at which particles switch from one state to the other
%spontaneously
r1=0.0225;

r2=0.0453;

% Define the range over which we will allow $r_3$ to vary to the left of
% the bifurcation point
r3_left=0:0.001:4*r1;
% Define the range over which we will allow $r_3$ to vary to the right of
% the bifurcation point
r3_right=4*r1:0.001:0.25;


%For these value of r_3 (either side of the bifurcation point) define the values of the steady state
z_plus=sqrt(1-4*r1./r3_right);

z_minus=-sqrt(1-4*r1./r3_right);

z0_stab=zeros(1,length(r3_left));

z0_unstab=zeros(1,length(r3_right));

%plot the steady states with their appropriate stability
plot(r3_right,z_plus,'linewidth',5)
hold on
plot(r3_right,z_minus,'linewidth',5)
plot(r3_left,z0_stab,'linewidth',5)
plot(r3_right,z0_unstab,'b--','linewidth',5)
%Define the xlabel
xlabel('r_3')

%save the figure
exportfig(gcf,...
        ['model_bifurcation_diagram.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['model_bifurcation_diagram.fig'],'fig');

