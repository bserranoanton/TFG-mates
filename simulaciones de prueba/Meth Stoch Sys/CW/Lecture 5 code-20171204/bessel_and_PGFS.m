%This code is desgined to plot some example Bessel functions and the
%stationary probability generating function with its derivative.
%for stochastic dimerisation decay and production reactions (4.1).
% This will reproduce figure 4.2 (a) and (b)
%By Kit Yates
%Created 27/09/15
%Last Modified 27/09/15

clear all
close all

%Define the rate constant of degradation
k1onu=0.005; %(sec^{-1})
k2nu=1; %(sec^{-1})

%%%%First reproduce the sample Bessel functions

% Define the range over which we will plot the Bessel Functions
z=0:0.01:4;

%Calculate the first 3 bessel functions
I0=besseli(0,z);
I1=besseli(1,z);
I2=besseli(2,z);

%plot the first three modified bessel functions of the first kind
plot(z,I0,'linewidth',5);
hold on
plot(z,I1,'r','linewidth',5);
plot(z,I2,'k','linewidth',5);
%Add a legend
legend('I_0(z)','I_1(z)','I_2(z)','location','NorthWest');
%Label the axes
xlabel('z')
ylabel('Modified Bessel functions of the first kind')

exportfig(gcf,...
            ['modified_bessel_functions.eps'],...
            'Format','eps2',...
            'Width','20',...
            'Color','cmyk',...
            'Resolution',300,...
            'FontMode','fixed',...
            'FontSize',21);
%             'LineWidth',6);
        %Save as a .fig as well
        saveas(gcf,['modified_bessel_functions.fig'],'fig');
        
        
%         Define the domain over which we will calculate the PGF
x=0:0.01:1;

% First calculate the normalisation constant C according to (4.19)
C=1/(sqrt(2)*besseli(1,2*sqrt(2*k2nu/k1onu)));
%The calculate the PGF according to (4.16)
G=C*sqrt(1+x).*besseli(1,2*sqrt(k2nu*(1+x)/k1onu));

%Calculate the derivative of the bessel function of the first kind
%according to the reccurance relation dI_v/dz=(v/z)*I_v(z)+I_{v+1}(z)
I1_prime=(1./(2*sqrt(k2nu*(1+x)/k1onu))).*besseli(1,2*sqrt(k2nu*(1+x)/k1onu))+besseli(2,2*sqrt(k2nu*(1+x)/k1onu));

%Finally calculate dG/dx
G_prime=(C./2*(sqrt(1+x))).*besseli(1,2*sqrt(k2nu*(1+x)/k1onu))+C*sqrt(k2nu/k1onu)*I1_prime;
  
%open a new figure
figure

%Plot G and its derivative on the same figure
[AX,H1,H2]=plotyy(x,G,x,G_prime);
  
%Add a legend on the top left corner
legend('G_s(x)','dG_s/dx(x)','location','NorthWest')
%Set the Ylimits of th eright hand axis
set(AX(2),'Ylim',[0,15])
%set the Yticks of the right hand axis
set(AX(2),'YTick',[0,5,10,15])
%Set the colors of the two axes
set(AX,{'ycolor'},{'b';'r'})
%Set the linewidth of the first plot
set(H1,'linewidth',5)
%Set the colro and with line width of the second plot
set(H2,'color','r','linestyle','--','linewidth',5)
%Set the ylabel for the left axis
set(get(AX(1),'Ylabel'),'String','G_s(x)')
%Set the ylabel for the right axis
set(get(AX(2),'Ylabel'),'String','dG_s/dx(x)')
%Det the xlabel
xlabel('x')

%Save the figure
exportfig(gcf,...
            ['generating_function_and_derivative.eps'],...
            'Format','eps2',...
            'Width','20',...
            'Color','cmyk',...
            'Resolution',300,...
            'FontMode','fixed',...
            'FontSize',21);
%             'LineWidth',6);
        %Save as a .fig as well
        saveas(gcf,['generating_function_and_derivative.fig'],'fig');