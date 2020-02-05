%This code is designed to plot the data generated in order to find the SPD
%of the minimal locusst model
%By Kit Yates
%Created 07/11/13
%Last modified 07/11/13

close all
clear all

%Find the name of the current directory and add it to the path
PRINTWORKINGDIRECTORY=pwd;
addpath(pwd);

%Decide whther or not we want to save the figures
SAVE_FIGURES=1;


%Define the name we will use to save the data
save_name='locusts_minimal_SPD';
directory_name='Figures/SPD';

%Define the total number of individuals
N=20;

%Augment the save name
save_name=[save_name,'N=',num2str(N)];
directory_name=[directory_name,'/N=',num2str(N)];

%Define the rate at which particles switch from one state to the other
%spontaneously
r1=0.0225;

%Augment the save name
save_name=[save_name,'r1=',num2str(r1)];
directory_name=[directory_name,'/r1=',num2str(r1)];

%Define the rates at which particles switch as a result of encounters
r2=0.0453;

%Augment the save name
save_name=[save_name,'r2=',num2str(r2)];
directory_name=[directory_name,'/r2=',num2str(r2)];

r3=0.1664;

%Augment the save name
save_name=[save_name,'r3=',num2str(r3)];
directory_name=[directory_name,'/r3=',num2str(r3)];

r4=0;

%Augment the save name
save_name=[save_name,'r4=',num2str(r4)];
directory_name=[directory_name,'/r4=',num2str(r4)];

%Add .mat to the save_name
save_name=[save_name,'.mat'];

%Load the relevant data
load(save_name);

%Determine the centres of those bins
bin_centres=2*(0:1:N)/N-1;

%Find the bin_spacing
bin_spacing=2/N;

%The histograms should integrate to 1 too
bar_vals=X1_rec/(sum(X1_rec)*bin_spacing);

%Symmetrise the data
bar_vals=(bar_vals+flipud(bar_vals))/2;

%plot the renormalised bar graph
bar(bin_centres,bar_vals,1);

%Define the spacing for the z points when plotting the SPD
z_spacing=0.001;

%define the spatial domain which we will use to plot the analytical SPD
z_values=-1:z_spacing:1;

%Define the values of the stationary probability distribution
% SPD=2*N*exp((N*z_values.^2)/2).*(4*r1+r3*(1-z_values.^2)).^(4*N*r1/r3-1);


%Define the groups of parameters that well be useful
a=8*r1*r4+4*r2^2+ 4*r2*(r3 + r4) + (r3+r4)^2;
b=4*r1*r4- 2*r2*(r3 + r4)-(r3+r4)^2;
c=r3+r4;
d=2*r2 + r3;

if r4~=0
    SPD=((d+r4*z_values.^2-sqrt(a)).^(N*(c*sqrt(a)+b)/((2*r4*sqrt(a))))).*((d+r4*z_values.^2+sqrt(a)).^(N*(c*sqrt(a)-b)/((2*r4*sqrt(a)))))./((r1+(2*r2+r4)*(1-z_values.^2)/4+r4*(1-z_values.^4)/8));
else
    SPD=((4*r1+(1-z_values.^2)*(2*r2+r3)).^(4*N*r1*(r2+r3)/(2*r2+r3)^2-1)).*exp((r3*N*(z_values.^2))/(2*(2*r2+r3)));
end

%Calculate the integral of the SPD
norm_factor=sum(SPD)*z_spacing;

%normalise the SPD
SPD_normed=SPD/norm_factor;

%Calculate also the SPD from the discrete-based model
SPD_discrete=zeros(N+1,1);

%Define the propensity functions as functions
T_plus = @(x) (1-x/N)*(r1+r2*x/N+r3*(x/N)^2+r4*(x/N)^3);
T_minus = @(x) (x/N)*(r1+r2*(1-x/N)+r3*(1-x/N)^2+r4*(1-x/N)^3);

%Initialise the prob of being in the first intervsl to be one (we will normalise later)
SPD_discrete(1)=1;

%Calculate the second value based on the boundary update equation
SPD_discrete(2)=SPD_discrete(1)*T_plus(0)/T_minus(1);

%Iteratively calculate the other value of 
for i=2:N
   SPD_discrete(i+1)=(SPD_discrete(i)*(T_plus(i-1)+T_minus(i-1))-SPD_discrete(i-1)*T_plus(i-2))/T_minus(i); 
end

% Finally normalise the probability distribution
SPD_discrete_normed=SPD_discrete/(sum(SPD_discrete)*bin_spacing);

%find the largest value of the deteministic SPD
max_SPD=max(SPD_normed);

%Find the largest value of the histogram
max_bar_vals=max(bar_vals);

%Find the largest value of the discrete probability distribution
max_SPD_discrete=max(SPD_discrete_normed)

%Find the height of the SPD plot
max_y=max([max_SPD,max_bar_vals,max_SPD_discrete]);


%plot also the analytically derived stationary distribution
hold on
plot(z_values,SPD_normed,'k','linewidth',4)
%Finally plot the normed discrete SPD as green stars
% plot(bin_centres,SPD_discrete_normed,'g*')
    xlabel('z','fontsize',30)
    ylabel('frequency','fontsize',30)
%Change the axis
axis([-1-bin_spacing/2 1+bin_spacing/2 0 max_y])
set(gca,'XTick',[-1 -0.5 0 0.5 1])
if SAVE_FIGURES
    
    %Make a directory with the relavant name
    mkdir(directory_name);
    
    set(gcf,'PaperUnits','centimeters');
    
    %Move to the correct directory
    cd(directory_name)
    
    exportfig(gcf,...
        [save_name,'.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',30,...
        'LineWidth',6);
    
    
    %Save the figure also as a jpg and a fig
    saveas(gcf,[save_name,'.fig'],'fig')
    saveas(gcf,[save_name,'.jpg'],'jpg')
    
    %cd back to the original directory
    cd(PRINTWORKINGDIRECTORY)
    
end