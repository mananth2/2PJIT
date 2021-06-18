%% %%%%%%%%%%%%%%%%% %%
%%%%% twoPJIT.m %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

% twoPJIT.m is the Main code that will be run by the user in order to
% obtain the dispersion plots and the perturbation quantities. The
% properties of the jet are entered in this script. This script calls
% Cylindrical_3D_Solution.m to descretize the system and solve for the
% eigenvalues. This script also calls perturbation.m in order to evaluate
% the perturbation quantities at the most unstable point.

%% %%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% List of variables %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ujet = liquid jet velocity (scalar)
% R = Radius of the jet (scalar)
% H = Gas domain boundary (scalar)
% rhol = density of the liquid jet (scalar)
% rhog = density of the gas (scalar)
% nul = kinematic viscosity of liquid (scalar)
% nug = kinematic viscosity of gas (scalar)
% delta_l = liquid shear layer thickness (scalar)
% delta_g = gas shear layer thickness (scalar)
% gamma = surface tension coefficient (scalar)
% m_theta = azimuthal wavenumber (0 for axisymmetric and 1 for asymmetric) (scalar)
% omegaArray = set of input non-dimensional frequency (=omega*R/Ujet) where
% omega is the dimensional frequency (vector)
% Nl = number of Gauss-Lobatto (G-L) points in liquid (scalar)
% Ng = number of G-L points in gas (scalar)
% Wavenumber = non-dimensional wavenumber obtained after calling
% Cylindrical_3D_solution (vector)
% GrowthRate = non-dimensional spatial growth rate obtained after calling
% Cylindrcal_3D_solution (vector)
% Cheb_coeff = eigenvector containing the chebyshev coefficients of
% perturbed velocities and pressure (2D-array)
% dominant_Grate = maximum non-dimensional growth rate obtained for the
% given set of input frequencies (scalar)
% dominant_wavenumber = non-dimensional wavenumber corresponding to the
% maximum growth rate (scalar)
% dominant_wavenumber = non-dimensional frequency corresponding to the
% maximum growth rate (scalar)
% ur_pert = perturbation velocity in radial direction obtained after
% calling perturbation.m (vector)
% utheta_pert = perturbation velocity in azimuthal direction obtained after
% calling perturbation.m (vector)
% uz_pert = perturbation velocity in axial direction after calling
% perturbation.m (vector)
% r = radial distance (vector)
% x, y and z = coordinates in 3D (vectors)
% X, Y and Z = coordinates in a 3D mesh (3D-arrays)
% interface = interface location for the given coordinates (scalar)
% radial = radial coordinates obtained from X and Y (scalar)
% alphaLiq = liquid volume fraction values in the 3D-mesh (3D-array)

clear all;
close all;
clc;

set(0,'defaultaxeslinewidth',3)
set(0,'defaultlinelinewidth',3)
set(0,'defaultlinemarkersize',10)
set(0,'defaultaxeslinewidth',3)
set(0,'defaultpatchlinewidth',3)
set(0,'defaultaxesfontsize',24)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Properties of the jet %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The properties of the jet are entered here
% The following values are default values which can be changed by the user

Ujet = 200;     
R = 45e-6;      
H = 5*R;        

rhol = 666.7;   
rhog = 50;      

nul = 6.947e-7; 
nug = 3.76e-7;  

delta_l = R/5;  
delta_g = R;    

gamma = 0.02;   


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Prompt for the mode of instability (Axisymmetric or Asymmetric) %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_theta = input(['Enter the mode of disturbance:'...
    '\n 0 - Axisymmetric mode'...
    '\n 1 - Asymmetric mode \n']);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Prompt for input frequency values %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% omegaArray here is the non-dimensional frequency given by:
% \tilde{omega}_R=omega_R*Ujet/R

% The prompt asks to enter values of the input frequencies which is entered
% as 'linspace' or 'logspace'. A default set 'linspace(0,10,40)' will be
% used. If the peak for the growth rate is observed well before
% end value of omegaArray, the end value of omegaArray can be reduced.
% Otherwise, if the peak for growth rate is not observed for a given
% omegaArray, the end value of omegaArray can be increased.

omegaArray=input(['\nInput an array of non-dimensional frequency values \n'...
    'Hit return to use a default value (linspace(0,10,40) will be used as default value) \n'...
    'Otherwise enter the values (usually linspace or logspace can be used) \n']);
if isempty(omegaArray)
    omegaArray = linspace(0,10,40);
end

omegaArray = omegaArray';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Prompt for the number of Gauss-Lobatto points %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In this section, number of Gauss-Lobatto (G-L) points in liquid and gas domains
% will be entered separately. Default number of points in liquid (Nl=40)
% and default number of points in gas (Ng=70) will be used. The code can
% be run for different values of Nl and Ng. It should be ensure that
% sufficient number of G-L points are used such that the dispersion curves
% converge as Nl and Ng are increased.

%% Number of G-L points in liquid

Nl=input(['\nEnter the number of Gauss-Lobatto points in liquid domain: \n'...
    'Hit return to use a default value (40 number of G-L points in liquid will be used as default value) \n'...
    'Otherwise enter the value \n']);
if isempty(Nl)
    Nl = 40;
end

%% Number of G-L points in gas
Ng=input(['Enter the number of Gauss-Lobatto points in gas domain: \n'...
    'Hit return to use a default value (70 number of G-L points in gas will be used as default value) \n'...
    'Otherwise enter the value \n']);
if isempty(Ng)
    Ng = 70;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Solving the instability by calling 'Cylindrical_3D_solution' function %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Wavenumber,GrowthRate,Cheb_coeff]=Cylindrical_3D_solution(Ujet,R,delta_l,delta_g,rhol,rhog,nul,nug,gamma,m_theta,omegaArray,H,Nl,Ng);
[M,I] = max(GrowthRate);  % M is the maxmimum growth rate value obtained for the given set of input frequencies. I is the corresponding index.

%% Most unstable gowth rate, wavenumber and the corresponding input frequency for the given set of input frequencies
dominant_Grate = GrowthRate(I,1);
dominant_wavenumber = Wavenumber(I,1);
dominant_frequency = omegaArray(I,1);


%% %%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Dispersion plots %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Disp = input(['\nSelect whether to plot non-dimensional or dimensional dispersion plot:'...
    '\n 0 - Non-dimensional'...
    '\n 1 - Dimensional \n']);
if (Disp == 0)  % Non-dimensional plots
    % Dispersion plot-1 (frequency, \omega_R vs growth rate, -k_I)
    f1=figure(1);
    movegui(f1,'north')
    plot(omegaArray,GrowthRate,'o--')
    xlabel('$Frequency, \hspace{0.2cm} \tilde{\omega}_R$','Interpreter','Latex')
    ylabel('$Growth \hspace{0.1cm} rate, \hspace{0.2cm} -\tilde{k}_I$','Interpreter','Latex')
    grid on
    box on
    hold on
        
    % Dispersion plot-2 (wavenumber, k_R vs growth rate, -k_I)
    f2=figure(2);
    movegui(f2,'north')
    plot(Wavenumber,GrowthRate,'o--')
    xlabel('$Wavenumber, \hspace{0.2cm} \tilde{k}_R$','Interpreter','Latex')
    ylabel('$Growth \hspace{0.1cm} rate, \hspace{0.2cm} -\tilde{k}_I$','Interpreter','Latex')
    grid on
    box on
    hold on
elseif (Disp == 1)  % Dimensional plots
    % Dispersion plot-1 (frequency, \omega_R in rad/s vs growth rate, -k_I in m^{-1})
    f1=figure(1);
    movegui(f1,'north')
    plot(omegaArray*Ujet/R,GrowthRate/R,'o--')
    xlabel('$Frequency, \hspace{0.2cm} \omega_R \hspace{0.2cm} (s^{-1})$','Interpreter','Latex')
    ylabel('$Growth \hspace{0.1cm} rate, \hspace{0.2cm} -k_I \hspace{0.2cm} (m^{-1})$','Interpreter','Latex')
    grid on
    box on
    hold on
        
    % Dispersion plot-2 (wavenumber, k_R in m^{-1} vs growth rate, -k_I in m^{-1})
    f2=figure(2);
    movegui(f2,'north')
    plot(Wavenumber/R,GrowthRate/R,'o--')
    xlabel('$Wavenumber, \hspace{0.2cm} k_R \hspace{0.2cm} (m^{-1})$','Interpreter','Latex')
    ylabel('$Growth \hspace{0.1cm} rate, \hspace{0.2cm} -k_I \hspace{0.2cm} (m^{-1})$','Interpreter','Latex')
    grid on
    box on
    hold on
else
    disp('Invalid input');
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['Frequency of the most unstable mode=', num2str(dominant_frequency*Ujet/R) ' rad/s']);
disp(['Wavenumber of the most unstable mode=', num2str(dominant_wavenumber/R) ' m^{-1}']);
disp(['Growth rate of the most unstable mode=', num2str(dominant_Grate/R) ' m^{-1}']);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Plots for profiles of perturbation quantties for highest growth rate mode %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'perturbation.m' is called to extract the perturbation quantities. These
% quantities are extracted for the omegaArray value for which the growth
% rate is the highest

[ur_pert,utheta_pert,uz_pert,p_pert,r]=perturbation(Ujet,R,delta_l,delta_g,rhol,rhog,nul,nug,gamma,m_theta,dominant_frequency,H,Nl,Ng);

% perturbation velocity in radial direction
f3=figure(3);
movegui(f3,'north')
plot(r,abs(ur_pert),'o--')
hold on
plot(R*ones(1,10),linspace(min(abs(ur_pert)),max(abs(ur_pert)),10),'k--')
ylabel('$|\hat{u}_r|, \hspace{0.2cm} (ms^{-1})$','Interpreter','Latex')
xlabel('$r \hspace{0.2cm} (m)$','Interpreter','Latex')
leg1 = legend(['$|\hat{u}_r|$'],['Unperturbed interface' char(10) '  location, $R$']);
set(leg1,'Interpreter','latex');
grid on
box on

% perturbation velocity in azimuthal direction
f4=figure(4);
movegui(f4,'north')
plot(r,abs(utheta_pert),'o--')
hold on
plot(R*ones(1,10),linspace(min(abs(utheta_pert)),max(abs(utheta_pert)),10),'k--')
ylabel('$|\hat{u}_{\theta}|, \hspace{0.2cm} (ms^{-1})$','Interpreter','Latex')
xlabel('$r \hspace{0.2cm} (m)$','Interpreter','Latex')
leg1 = legend(['$|\hat{u}_{\theta}|$'],['Unperturbed interface' char(10) '  location, $R$']);
set(leg1,'Interpreter','latex');
grid on
box on

% perturbation velocity in axial direction
f5=figure(5);
movegui(f5,'north')
plot(r,abs(uz_pert),'o--')
hold on
plot(R*ones(1,10),linspace(min(abs(uz_pert)),max(abs(uz_pert)),10),'k--')
ylabel('$|\hat{u}_z|, \hspace{0.2cm} (ms^{-1})$','Interpreter','Latex')
xlabel('$r \hspace{0.2cm} (m)$','Interpreter','Latex')
leg1 = legend(['$|\hat{u}_z|$'],['Unperturbed interface' char(10) '  location, $R$']);
set(leg1,'Interpreter','latex');
grid on
box on

% perturbation pressure
f6=figure(6);
movegui(f6,'north')
plot(r,abs(p_pert),'o--')
hold on
plot(R*ones(1,10),linspace(min(abs(p_pert)),max(abs(p_pert)),10),'k--')
ylabel('$|\hat{p}|, \hspace{0.2cm} (Nm^{-2})$','Interpreter','Latex')
xlabel('$r \hspace{0.2cm} (m)$','Interpreter','Latex')
leg1 = legend(['$|\hat{p}|$'],['Unperturbed interface' char(10) '  location, $R$']);
set(leg1,'Interpreter','latex');
grid on
box on


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Visual representation of the dominant mode %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to plot the cylindrical jet with the dominant mode of
% perturbation obtained for the given set of input frequencies.

x=linspace(-2*R,2*R,100);
y=linspace(-2*R,2*R,100);
z=linspace(0,20*R,200);
[Z,X,Y]=meshgrid(z,x,y);
for i=1:size(Z,1)
for j=1:size(Z,2)
for k=1:size(Z,3)
    if (X(i,j,k)>=0)
        interface(i,j,k)=R+0.05*R*exp(dominant_Grate*Z(i,j,k)/R)*cos(dominant_wavenumber*Z(i,j,k)/R+m_theta*atan((Y(i,j,k)/X(i,j,k))));
    else
        interface(i,j,k)=R+0.05*R*exp(dominant_Grate*Z(i,j,k)/R)*cos(dominant_wavenumber*Z(i,j,k)/R+m_theta*(atan((Y(i,j,k)/X(i,j,k)))+pi));
    end
radial(i,j,k)=sqrt((X(i,j,k))^2+(Y(i,j,k))^2);
if((radial(i,j,k))>(interface(i,j,k)))
alphaLiq(i,j,k)=0;
else
alphaLiq(i,j,k)=1;
end
end
end
end
alphaLiq_smooth = smooth3(alphaLiq);

figure(7)
p = patch(isosurface(Z,X,Y,alphaLiq_smooth,0.5));
isonormals(Z,X,Y,alphaLiq_smooth,p)
p.FaceColor = '[0 0.5 1]';
p.EdgeColor = 'none';
daspect([1 1 1])
view(37.5,30); 
axis tight 
camlight
lighting gouraud
xlabel('$z$','Interpreter','Latex')
ylabel('$x$','Interpreter','Latex')
zlabel('$y$','Interpreter','Latex')