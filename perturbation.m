%% %%%%%%%%%%%%%%%%%%%%%% %%
%%%%% perturbation.m %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perturbation.m evaluates the perturbation velocities and pressure at the
% most unstable mode by calling 'Cylindrical_3D_Solution.m'

%% %%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% List of variables %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ujet = liquid jet velocity (scalar) (m/s)
% R = Radius of the jet (scalar) (m)
% H = Gas domain boundary (scalar) (m)
% rho_L = density of the liquid jet (scalar) (kg/m^3)
% rho_G = density of the gas (scalar) (kg/m^3)
% nu_L = kinematic viscosity of liquid (scalar) (m^2/s)
% nu_G = kinematic viscosity of gas (scalar) (m^2/s)
% delta_L = liquid shear layer thickness (scalar) (m)
% delta_G = gas shear layer thickness (scalar) (m)
% gamma = surface tension coefficient (scalar) (N/m)
% m_theta = azimuthal wavenumber (0 for axisymmetric and 1 for asymmetric) (scalar)
% omega = input non-dimensional frequency (scalar)
% N_L = number of Gauss-Lobatto (G-L) points in liquid (scalar)
% N_G = number of G-L points in gas (scalar)

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
% vGamma = non-dimensional radial velocity at the interface used for scaling (scalar) 
% ur_pert = perturbation velocity in radial direction (vector) (m/s)
% utheta_pert = perturbation velocity in azimuthal direction (vector) (m/s)
% uz_pert = perturbation velocity in axial direction (vector) (m/s)
% r_full = radial distance (vector) (m)

function[ur_pert,utheta_pert,uz_pert,p_pert,r_full]=perturbation(Ujet,R,delta_L,delta_G,rho_L,rho_G,nu_L,nu_G,gamma,m_theta,omega,H,N_L,N_G)

[Wavenumber,GrowthRate,Cheb_coeff]=Cylindrical_3D_solution(Ujet,R,delta_L,delta_G,rho_L,rho_G,nu_L,nu_G,gamma,m_theta,omega,H,N_L,N_G);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Defining Chebyshev matrices %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Liquid domain
% Chebyshev matrices are defined for internal liquid G-L points excluding
% liquid centerline (r=0 or \tilde{r}_L=-1} and
% interface on liquid side (r=R or \tilde{r}_L=+1)
% These are used for all liquid governing equations

for j=1:N_L-1
    rL(j,1) = -cos((pi/((N_L)))*j);    % internal liquid G-L points
end

vec_liq = acos(rL);

% The Chebyshev matrices are defined below
D0_liq=[];
vec=(1:1:N_L-1)';
for j=0:1:N_L
    D0_liq=[D0_liq cos(j*vec_liq)];
end

%% Gas domain
% Chebyshev matrices are defined for internal gas G-L points excluding
% gas boundary (r=L or \tilde{r}_G=+1} and
% interface on gas side (r=R or \tilde{r}_G=-1)
% These are used for all gas governing equations

for j=1:N_G-1
    rG(j,1) = -cos((pi/((N_G)))*j);  % internal gas G-L points
end
vec_gas = acos(rG);

% The Chebyshev matrices are defined below 
D0_gas=[];
for j=0:1:N_G
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
for j=0:1:N_L
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
for j=0:1:N_G
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
for j=0:1:N_L
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
for j=0:1:N_G
    D0_gasInterface=[D0_gasInterface cos(j*vec_gasInterface)];
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Determining perturbation quantities using Chebyshev series %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vGamma = D0_liqInterface*Cheb_coeff(1:N_L+1,1);
epsilon = 0.01/abs(vGamma);
ur_pert = Ujet*[D0_liq*Cheb_coeff(1:N_L+1,1); D0_gas*Cheb_coeff(4*N_L+4:4*N_L+3+N_G+1,1)]*epsilon;
utheta_pert = Ujet*[D0_liq*Cheb_coeff(N_L+2:2*N_L+2,1); D0_gas*Cheb_coeff(4*N_L+3+N_G+2:4*N_L+3+2*N_G+2,1)]*epsilon;
uz_pert = Ujet*[D0_liq*Cheb_coeff(2*N_L+3:3*N_L+3,1); D0_gas*Cheb_coeff(4*N_L+3+2*N_G+3:4*N_L+3+3*N_G+3,1)]*epsilon;
p_pert = rho_L*Ujet^2*[D0_liq(:,1:end-1)*Cheb_coeff(3*N_L+4:4*N_L+3,1); D0_gas(:,1:end-1)*Cheb_coeff(4*N_L+3+3*N_G+4:4*N_L+3+4*N_G+3,1)]*epsilon;
r_full = [R/2*(rL+1); R/2*(rG*(H/R-1)+H/R+1)];

end