%% %%%%%%%%%%%%%%%%%%%%%% %%
%%%%% perturbation.m %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perturbation.m evaluates the perturbation velocities and pressure at the
% most unstable mode by calling 'Cylindrical_3D_Solution.m'

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
% omega = input non-dimensional frequency (scalar)
% Nl = number of Gauss-Lobatto (G-L) points in liquid (scalar)
% Ng = number of G-L points in gas (scalar)

% rL = transformed radial positions of internal liquid G-L points (vector)
% rG = transformed radial positions of internal gas G-L points (vector)
% rL_liqCenter = transformed radial position for liquid centerline point (scalar)
% rG_gasBoundary = transformed radial position for gas boundary point (scalar)
% rL_interface = transformed radial position for liquid interface point (scalar)
% rG_interface = tranformed radial position for gas interface point (scalar)

% D0_liq = Chebyshev matrix for internal liquid G-L points (2D-array)
% D0_gas = Chebyshev matrix for internal gas G-L points (2D-array)
% D0_liqCenter = Chebyshev matrix for liquid centerline point (row-matrix)
% D0_gasBoundary = Chebyshev matrix for gas boundary point (row-matrix)
% D0_liqInterface = Chebyshev matrix for liquid interface point (row-matrix)
% D0_gasInterfacee = Chebyshev matrix for gas interface point (row-matrix)

% Cheb_coeff = eigenvector containing the chebyshev coefficients of
% perturbed velocities and pressure (vector)
% Wavenumber = non-dimensional wavenumber obtained (scalar)
% Growth rate = non-dimensional growth rate obtained (scalar)

% epsilon = scaling down factor for the eigenvector (scalar)
% ur_pert = perturbation velocity in radial direction (vector)
% utheta_pert = perturbation velocity in azimuthal direction (vector)
% uz_pert = perturbation velocity in axial direction (vector)
% r_full = radial distance (vector)

function[ur_pert,utheta_pert,uz_pert,p_pert,r_full]=perturbation(Ujet,R,delta_l,delta_g,rhol,rhog,nul,nug,gamma,m_theta,omega,H,Nl,Ng)

[Wavenumber,GrowthRate,Cheb_coeff]=Cylindrical_3D_solution(Ujet,R,delta_l,delta_g,rhol,rhog,nul,nug,gamma,m_theta,omega,H,Nl,Ng);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Defining Chebyshev matrices %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Liquid domain
% Chebyshev matrices are defined for internal liquid G-L points excluding
% liquid centerline (r=0 or \tilde{r}_L=-1} and
% interface on liquid side (r=R or \tilde{r}_L=+1)
% These are used for all liquid governing equations

for j=1:Nl-1
    rL(j,1) = -cos((pi/((Nl)))*j);    % internal liquid G-L points
end

vec_liq = acos(rL);

% The Chebyshev matrices are defined below
D0_liq=[];
vec=(1:1:Nl-1)';
for j=0:1:Nl
    D0_liq=[D0_liq cos(j*vec_liq)];
end

%% Gas domain
% Chebyshev matrices are defined for internal gas G-L points excluding
% gas boundary (r=L or \tilde{r}_G=+1} and
% interface on gas side (r=R or \tilde{r}_G=-1)
% These are used for all gas governing equations

for j=1:Ng-1
    rG(j,1) = -cos((pi/((Ng)))*j);  % internal gas G-L points
end
vec_gas = acos(rG);

% The Chebyshev matrices are defined below 
D0_gas=[];
for j=0:1:Ng
    D0_gas=[D0_gas cos(j*vec_gas)];
end

%% Liquid centerline
% Chebyshev matrices are defined for liquid for
% liquid centerline (r=0 or \tilde{r}_L=-1} 
% These are used for all centerline conditions

vec_liqCenter = pi;   
rL_liqCenter = cos(vec_liqCenter); % liquid centerline location
    
% The Chebyshev matrices are defined below
D0_liqCenter=[];
for j=0:1:Nl
    D0_liqCenter=[D0_liqCenter cos(j*vec_liqCenter)];
end
    
%% Gas boundary
% Chebyshev matrices are defined for 
% gas boundary (r=L or \tilde{r}_G=+1}
% These are used for all gas boundary conditions

vec_gasBoundary = 0;   
rG_gasBoundary = cos(vec_gasBoundary); % gas boundary location
    
% The Chebyshev matrices are defined below 
D0_gasBoundary=[];
for j=0:1:Ng
    D0_gasBoundary=[D0_gasBoundary cos(j*vec_gasBoundary)];
end
    
%% Liquid Interfacial location      
% Chebyshev matrices are defined for 
% liquid interface (r=R or \tilde{r}_L=+1}
% These used for all interface conditions

vec_liqInterface = 0;   
rL_interface = cos(vec_liqInterface);   % liquid interface location
    
% The Chebyshev matrices are defined below 
D0_liqInterface=[];
for j=0:1:Nl
    D0_liqInterface=[D0_liqInterface cos(j*vec_liqInterface)];
end

%% Gas Interfacial location     
% Chebyshev matrices are defined for 
% gas interface (r=R or \tilde{r}_G=-1}
% These are used for all interface conditions
    
vec_gasInterface = pi;
rG_gasInterface = cos(vec_gasInterface);   % gas interface location
 
% The Chebyshev matrices are defined below
D0_gasInterface=[];
for j=0:1:Ng
    D0_gasInterface=[D0_gasInterface cos(j*vec_gasInterface)];
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Determining perturbation quantities using Chebyshev series %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon = 1e-2;
ur_pert = Ujet*[D0_liq*Cheb_coeff(1:Nl+1,1); D0_gas*Cheb_coeff(4*Nl+4:4*Nl+3+Ng+1,1)]*epsilon;
utheta_pert = Ujet*[D0_liq*Cheb_coeff(Nl+2:2*Nl+2,1); D0_gas*Cheb_coeff(4*Nl+3+Ng+2:4*Nl+3+2*Ng+2,1)]*epsilon;
uz_pert = Ujet*[D0_liq*Cheb_coeff(2*Nl+3:3*Nl+3,1); D0_gas*Cheb_coeff(4*Nl+3+2*Ng+3:4*Nl+3+3*Ng+3,1)]*epsilon;
p_pert = rhol*Ujet^2*[D0_liq(:,1:end-1)*Cheb_coeff(3*Nl+4:4*Nl+3,1); D0_gas(:,1:end-1)*Cheb_coeff(4*Nl+3+3*Ng+4:4*Nl+3+4*Ng+3,1)]*epsilon;
r_full = [R/2*(rL+1); R/2*(rG*(H/R-1)+H/R+1)];

end