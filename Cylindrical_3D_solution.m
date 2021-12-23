%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Cylindrical_3D_Solution.m %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cylindrical_3D_Solution.m is used to discretize the system using 
% Chebyshev Spetral method and solving the final discretized eigenvalue 
% problem.
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
% omegaArray = set of input non-dimensional frequency (=omega*R/Ujet) where
% omega (rad/s) is the dimensional frequency (vector)
% N_L = number of Gauss-Lobatto (G-L) points in liquid (scalar)
% N_G = number of G-L points in gas (scalar)
% b = Gas domain length in radial direction (=H-R) (scalar) (m)
% lr = Ratio of total domain length in radial direction to the radius of
% the jet (scalar)
% Dr = density ratio (scalar)
% mul = dynamic viscosity of liquid (scalar) (N-s/m^2)
% mug = dynamic viscosity of gas (scalar) (N-s/m^2)
% Vr = viscosity ratio (scalar)
% Ugstar and Ulstar = velocities used in the evaluation of error function
% profile for base velocity (scalars) (m/s)
% Rel = Reynolds number using Ujet, rhol, R and mul (scalar)
% Wel = Weber number using rhol, Ujet, R and gamma (scalar)
% rL = transformed radial positions of internal liquid G-L points (vector) 
% rG = transformed radial positions of internal gas G-L points (vector)
% rL_liqCenter = transformed radial position for liquid centerline point (scalar)
% rG_gasBoundary = transformed radial position for gas boundary point (scalar)
% rL_interface = transformed radial position for liquid interface point (scalar)
% rG_interface = tranformed radial position for gas interface point (scalar)

% D0_liq, D1_liq and D2_liq = Chebyshev matrix, first derivative and second
% derivative of Chebyshev matrix for internal liquid G-L points (2D-arrays)
% D0_gas, D1_gas and D2_gas = Chebyshev matrix, first derivative and second
% derivative of Chebyshev matrix for internal gas G-L points (2D-arrays)
% D0_liqCenter, D1_liqCenter and D2_liqCenter = Chebyshev matrix, first derivative and second
% derivative of Chebyshev matrix for liquid centerline point (row matrices)
% D0_gasBoundary, D1_gasBoundary and D2_gasBoundary = Chebyshev matrix, first derivative and second
% derivative of Chebyshev matrix for gas boundary point (row matrices)
% D0_liqInterface, D1_liqInterface and D2_liqInterface = Chebyshev matrix, first derivative and second
% derivative of Chebyshev matrix for liquid interface point (row matrices)
% D0_gasInterface, D1_gasInterface and D2_gasInterface = Chebyshev matrix, first derivative and second
% derivative of Chebyshev matrix for gas interface point (row matrices)

% Ul = non-dimensional velocity at liquid internal G-L points (vector)
% Ulp = first derivative of non-dimensional velocity at liquid internal G-L
% internal points with respect to transformed radial distance rL (vector)
% Ulpp = second derivative of non-dimensional velocity at liquid internal G-L
% internal points with respect to transformed radial distance rL (vector)

% Ug = non-dimensional velocity at gas internal G-L points (vector)
% Ugp = first derivative of non-dimensional velocity at gas internal G-L
% internal points with respect to transformed radial distance rG (vector)
% Ugpp = second derivative of non-dimensional velocity at gas internal G-L
% internal points with respect to transformed radial distance rG (vector)

% C0, C1, C2 = corresponding matrices in the non-linear equation of the
% form k^2*C2.a+k*C1.a+C0.a=0 (2D-arrays)
% A, B = matrices in the linearized eigenvalue problem of the form A.f=k*B.f
% where f=[ka;a]' (2D-array)
% k = set of eigenvalues obtained by solving the eigenvalue problem (vector)
% F = set of eigenvectors obtained by solving the eigenvalue problem (vector)
% kFinal = Eigenvalue with maximum growth rate selected from the
% non-spurious eigenvalues (vector)
% Cheb_coeff = eigenvector corresponding to kFinal (eigenvalue of interest)
% (This contains chebyshev coefficients for the perturbed velocities and
% pressure) (2D-array)
% Wavenumber = non-dimensional wavenumber obtained (real part of kFinal) (vector)
% Growth rate = non-dimensional growth rate obtained (negative of the the
% imaginary part of kFinal) (vector)
% Error = Residual error for the linearized eigenvalue
% problem (||A.f-k*B.f||) (vector)



function[Wavenumber,GrowthRate,Cheb_coeff]=Cylindrical_3D_solution(Ujet,R,delta_L,delta_G,rho_L,rho_G,nu_L,nu_G,gamma,m_theta,omegaArray,H,N_L,N_G)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Properties of the jet %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = H-R;    

Dr = rho_G/rho_L;     
lr = H/R;       

mu_L = rho_L*nu_L;  
mu_G = rho_G*nu_G; 
Vr = mu_G/mu_L;   

Ugstar = Ujet/(1+Vr*delta_L/delta_G);   
Ulstar = Ugstar*Vr*delta_L/delta_G;     

Rel = Ujet*R/(nu_L);     
Wel = rho_L*Ujet^2*R/(gamma);  

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
lv_liq=length(vec);
D1_liq=[zeros(lv_liq,1) D0_liq(:,1) 4*D0_liq(:,2)];
D2_liq=[zeros(lv_liq,1) zeros(lv_liq,1) 4*D0_liq(:,1)];
for j=3:N_L
    D1_liq=[D1_liq 2*j*D0_liq(:,j)+j*D1_liq(:,j-1)/(j-2)];
    D2_liq=[D2_liq 2*j*D1_liq(:,j)+j*D2_liq(:,j-1)/(j-2)];
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
lv_gas=length(vec_gas);
D1_gas=[zeros(lv_gas,1) D0_gas(:,1) 4*D0_gas(:,2)];
D2_gas=[zeros(lv_gas,1) zeros(lv_gas,1) 4*D0_gas(:,1)];
for j=3:N_G
    D1_gas=[D1_gas 2*j*D0_gas(:,j)+j*D1_gas(:,j-1)/(j-2)];
    D2_gas=[D2_gas 2*j*D1_gas(:,j)+j*D2_gas(:,j-1)/(j-2)];
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
lv_liqCenter=length(vec_liqCenter);
D1_liqCenter=[zeros(lv_liqCenter,1) D0_liqCenter(:,1) 4*D0_liqCenter(:,2)];
D2_liqCenter=[zeros(lv_liqCenter,1) zeros(lv_liqCenter,1) 4*D0_liqCenter(:,1)];
for j=3:N_L
    D1_liqCenter=[D1_liqCenter 2*j*D0_liqCenter(:,j)+j*D1_liqCenter(:,j-1)/(j-2)];
    D2_liqCenter=[D2_liqCenter 2*j*D1_liqCenter(:,j)+j*D2_liqCenter(:,j-1)/(j-2)];
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
lv_gasBoundary=length(vec_gasBoundary);
D1_gasBoundary=[zeros(lv_gasBoundary,1) D0_gasBoundary(:,1) 4*D0_gasBoundary(:,2)];
D2_gasBoundary=[zeros(lv_gasBoundary,1) zeros(lv_gasBoundary,1) 4*D0_gasBoundary(:,1)];
for j=3:N_G
    D1_gasBoundary=[D1_gasBoundary 2*j*D0_gasBoundary(:,j)+j*D1_gasBoundary(:,j-1)/(j-2)];
    D2_gasBoundary=[D2_gasBoundary 2*j*D1_gasBoundary(:,j)+j*D2_gasBoundary(:,j-1)/(j-2)];
end
    
%% Liquid InterfacIAL LOCATION  
% Chebyshev matrices are defined for 
% liquid interface (r=R or \tilde{r}_L=+1}
% These are used for all interface conditions

vec_liqInterface = 0;   
rL_interface = cos(vec_liqInterface); % liquid interface location

% The Chebyshev matrices are defined below
D0_liqInterface=[];
for j=0:1:N_L
    D0_liqInterface=[D0_liqInterface cos(j*vec_liqInterface)];
end
lv_liqInterface=length(vec_liqInterface);
D1_liqInterface=[zeros(lv_liqInterface,1) D0_liqInterface(:,1) 4*D0_liqInterface(:,2)];
D2_liqInterface=[zeros(lv_liqInterface,1) zeros(lv_liqInterface,1) 4*D0_liqInterface(:,1)];
for j=3:N_L
    D1_liqInterface=[D1_liqInterface 2*j*D0_liqInterface(:,j)+j*D1_liqInterface(:,j-1)/(j-2)];
    D2_liqInterface=[D2_liqInterface 2*j*D1_liqInterface(:,j)+j*D2_liqInterface(:,j-1)/(j-2)];
end

%% Gas InterfacIAL LOCATION      
% Chebyshev matrices are defined for 
% gas interface (r=R or \tilde{r}_G=-1}
% These are used for all interface conditions
    
vec_gasInterface = pi;
rG_gasInterface = cos(vec_gasInterface); %gas interface location  

% The Chebyshev matrices are defined below
D0_gasInterface=[];
for j=0:1:N_G
    D0_gasInterface=[D0_gasInterface cos(j*vec_gasInterface)];
end
lv_gasInterface=length(vec_gasInterface);
D1_gasInterface=[zeros(lv_gasInterface,1) D0_gasInterface(:,1) 4*D0_gasInterface(:,2)];
D2_gasInterface=[zeros(lv_gasInterface,1) zeros(lv_gasInterface,1) 4*D0_gasInterface(:,1)];
for j=3:N_G
    D1_gasInterface=[D1_gasInterface 2*j*D0_gasInterface(:,j)+j*D1_gasInterface(:,j-1)/(j-2)];
    D2_gasInterface=[D2_gasInterface 2*j*D1_gasInterface(:,j)+j*D2_gasInterface(:,j-1)/(j-2)];
end

%% %%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Base velocities %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Liquid

% Base velocity (non-dimensionalized using Ujet)

Ul = -Ulstar*erf((rL-1)*R/(2*delta_L))+Ugstar;
Ul = Ul/Ujet;

Ul_int = -Ulstar*erf((1-1)*R/(2*delta_L))+Ugstar;
Ul_int = Ul_int/Ujet;     % interface value

% First derivative of velocity with respect to r_L (non-dimensionalized
% using Ujet)

Ulp = -Ulstar*(2/pi^(1/2))*(R/(2*delta_L))*exp(-((rL-1)*R/(2*delta_L)).^2);
Ulp = Ulp/Ujet;

Ulp_int = -Ulstar*(2/pi^(1/2))*(R/(2*delta_L))*exp(-((1-1)*R/(2*delta_L)).^2);
Ulp_int = Ulp_int/Ujet;   % interface value

% Second derivative of velocity with respect to r_L (non-dimensionalized
% using Ujet)

Ulpp = Ulstar*(R^2/(delta_L^2*pi^(1/2)))*((rL-1)*R/(2*delta_L)).*exp(-((rL-1)*R/(2*delta_L)).^2);
Ulpp = Ulpp/Ujet;

Ulpp_int = Ulstar*(R^2/(delta_L^2*pi^(1/2)))*((1-1)*R/(2*delta_L)).*exp(-((1-1)*R/(2*delta_L)).^2);
Ulpp_int = Ulpp_int/Ujet; % interface value

%% Gas

% Base velocity (non-dimensionalized
% using Ujet)

Ug = -Ugstar*erf((rG+1)*b/(2*delta_G))+Ugstar;
Ug = Ug/Ujet;

Ug_int = -Ugstar*erf((-1+1)*b/(2*delta_G))+Ugstar;
Ug_int = Ug_int/Ujet;     % interface value

% First derivative of velocity with respect to r_G (non-dimensionalized
% using Ujet)

Ugp = -Ugstar*(2/pi^(1/2))*(b/(2*delta_G))*exp(-((rG+1)*b/(2*delta_G)).^2);
Ugp = Ugp/Ujet;

Ugp_int = -Ugstar*(2/pi^(1/2))*(b/(2*delta_G))*exp(-((-1+1)*b/(2*delta_G)).^2);
Ugp_int = Ugp_int/Ujet;   % interface value

% Second derivative of velocity with respect to r_G (non-dimensionalized
% using Ujet)

Ugpp = Ugstar*(b^2/(delta_G^2*pi^(1/2)))*((rG+1)*b/(2*delta_G)).*exp(-((rG+1)*b/(2*delta_G)).^2);
Ugpp = Ugpp/Ujet;

Ugpp_int = Ugstar*(b^2/(delta_G^2*pi^(1/2)))*((-1+1)*b/(2*delta_G)).*exp(-((-1+1)*b/(2*delta_G)).^2);
Ugpp_int = Ugpp_int/Ujet; % interface value


Uint = Ul_int;  % There is continuity of base velocity at interface (Uint=Ul_int=Ug_int)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% The following section is solving the system of governing %%%%%%%
%%%%% equations along with boundary and interfacial conditions %%%%%%%
%%%%% for each frequency (omega) values. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% System of linear equations in Chebyshev space
for i=1:length(omegaArray)
    omega = omegaArray(i,1);    % This is the non-dimensional frequency
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Construction of matrices (C0,C1,C2) in the eigenvalue %%%%%%%
    %%%%% problem of the form: k^2*C2.a+k*C1.a+C0.a=0 %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Governing equations for liquid phase %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Continuity
    
    % C0 matrix
    C0_liq_11 = 2*D1_liq + (2./(rL+1)).*D0_liq; % radial velocity part
    C0_liq_12 = (2i*m_theta./(rL+1)).*D0_liq;   % azimuthal velocity part
    C0_liq_13 = zeros(N_L-1,N_L+1);               % axial velocity part
    C0_liq_14 = zeros(N_L-1,N_L+1);               % pressure part
    
    % C1 matrix
    C1_liq_11 = zeros(N_L-1,N_L+1);               % radial velocity part
    C1_liq_12 = zeros(N_L-1,N_L+1);               % azimuthal velocity part
    C1_liq_13 = 1i*D0_liq;                      % axial velocity part
    C1_liq_14 = zeros(N_L-1,N_L+1);               % pressure part
    
    % C2 matrix
    C2_liq_11 = zeros(N_L-1,N_L+1);               % radial velocity part
    C2_liq_12 = zeros(N_L-1,N_L+1);               % azimuthal velocity part
    C2_liq_13 = zeros(N_L-1,N_L+1);               % axial velocity part
    C2_liq_14 = zeros(N_L-1,N_L+1);               % pressure part
    
    %% radial momentum
    
    % C0 matrix
    C0_liq_21 = (4/Rel)*D2_liq + (4./(Rel*(rL+1))).*D1_liq +(-((4*m_theta^2)./(Rel*(rL+1).^2)) -(4./(Rel*(rL+1).^2)) + 1i*omega).*D0_liq;   % radial velocity part
    C0_liq_22 = -(1/Rel)*(1i*8*m_theta./((rL+1).^2)).*D0_liq;   % azimuthal velocity part
    C0_liq_23 = zeros(N_L-1,N_L+1);                               % axial velocity part
    C0_liq_24 = -2*D1_liq;                                      % pressure part
    
    % C1 matrix
    C1_liq_21 = -1i*Ul.*D0_liq;     % radial velocity part
    C1_liq_22 = zeros(N_L-1,N_L+1);   % azimuthal velocity part
    C1_liq_23 = zeros(N_L-1,N_L+1);   % axial velocity part
    C1_liq_24 = zeros(N_L-1,N_L+1);   % pressure part
    
    % C2 matrix
    C2_liq_21 = -(1/Rel)*D0_liq;    % radial velocity part
    C2_liq_22 = zeros(N_L-1,N_L+1);   % azimuthal velocity part
    C2_liq_23 = zeros(N_L-1,N_L+1);   % axial velocity part
    C2_liq_24 = zeros(N_L-1,N_L+1);   % pressure part
    
    %% azimuthal momentum
    
    % C0 matrix
    C0_liq_31 = (1/Rel)*(8*1i*m_theta./((rL+1).^2)).*D0_liq;    % radial velocity part
    C0_liq_32 = (4/Rel)*D2_liq + (4./(Rel*(rL+1))).*D1_liq +(-((4*m_theta^2)./(Rel*(rL+1).^2)) -(4./(Rel*(rL+1).^2)) + 1i*omega).*D0_liq;   % azimuthal velocity part
    C0_liq_33 = zeros(N_L-1,N_L+1);                   % axial velocity part
    C0_liq_34 = -(2*1i*m_theta./((rL+1))).*D0_liq;  % pressure part
    
    % C1 matrix
    C1_liq_31 = zeros(N_L-1,N_L+1);   % radial velocity part
    C1_liq_32= -1i*Ul.*D0_liq;      % azimuthal velocity part
    C1_liq_33 = zeros(N_L-1,N_L+1);   % axial velocity part
    C1_liq_34 = zeros(N_L-1,N_L+1);   % pressure part
    
    % C2 matrix
    C2_liq_31 = zeros(N_L-1,N_L+1);   % radial velocity part
    C2_liq_32= -(1/Rel)*D0_liq;     % azimuthal velocity part
    C2_liq_33 = zeros(N_L-1,N_L+1);   % axial velocity part
    C2_liq_34 = zeros(N_L-1,N_L+1);   % pressure part
    
    %% axial momentum
    
    % C0 matrix
    C0_liq_41 = -2*Ulp.*D0_liq;     % radial velocity part
    C0_liq_42 = zeros(N_L-1,N_L+1);   % azimuthal velocity part
    C0_liq_43 = (4/Rel)*D2_liq + (4./(Rel*(rL+1))).*D1_liq +(-((4*m_theta^2)./(Rel*(rL+1).^2)) + 1i*omega).*D0_liq; % axial velocity part
    C0_liq_44 = zeros(N_L-1,N_L+1);   % pressure part
    
    % C1 matrix
    C1_liq_41 = zeros(N_L-1,N_L+1);   % radial velocity part
    C1_liq_42 = zeros(N_L-1,N_L+1);   % azimuthal velocity part
    C1_liq_43 = -1i*Ul.*D0_liq;     % axial velocity part
    C1_liq_44 = -1i*D0_liq;         % pressure part
    
    % C2 matrix
    C2_liq_41 = zeros(N_L-1,N_L+1);   % radial velocity part
    C2_liq_42 = zeros(N_L-1,N_L+1);   % azimuthal velocity part
    C2_liq_43 = -(1/Rel)*D0_liq;    % axial velocity part
    C2_liq_44 = zeros(N_L-1,N_L+1);   % pressure part
    
    %% Construction of the C0,C1 and C1 matrices for liquid 
    
    % zero rows at the end is to fill the liquid centerline and interface
    % conditions. 
    C0_liq = [C0_liq_21 C0_liq_22 C0_liq_23 C0_liq_24(:,1:end-1);...
        C0_liq_31 C0_liq_32 C0_liq_33 C0_liq_34(:,1:end-1);...
        C0_liq_41 C0_liq_42 C0_liq_43 C0_liq_44(:,1:end-1);...
        C0_liq_11 C0_liq_12 C0_liq_13 C0_liq_14(:,1:end-1); zeros(7,4*N_L+3)];
    
    C1_liq = [C1_liq_21 C1_liq_22 C1_liq_23 C1_liq_24(:,1:end-1);...
        C1_liq_31 C1_liq_32 C1_liq_33 C1_liq_34(:,1:end-1);
        C1_liq_41 C1_liq_42 C1_liq_43 C1_liq_44(:,1:end-1);...
        C1_liq_11 C1_liq_12 C1_liq_13 C1_liq_14(:,1:end-1); zeros(7,4*N_L+3)];
    
    C2_liq = [C2_liq_21 C2_liq_22 C2_liq_23 C2_liq_24(:,1:end-1);...
        C2_liq_31 C2_liq_32 C2_liq_33 C2_liq_34(:,1:end-1);...
        C2_liq_41 C2_liq_42 C2_liq_43 C2_liq_44(:,1:end-1);...
        C2_liq_11 C2_liq_12 C2_liq_13 C2_liq_14(:,1:end-1); zeros(7,4*N_L+3)];
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Governing equations for gas phase %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %% Continuity
    
    % C0 matrix
    C0_gas_11 = (2/(lr-1))*D1_gas + (2./(rG*(lr-1)+lr+1)).*D0_gas;	% radial velocity part
    C0_gas_12 = (2i*m_theta./(rG*(lr-1)+lr+1)).*D0_gas; % azimuthal velocity part
    C0_gas_13 = zeros(N_G-1,N_G+1);   % axial velocity part
    C0_gas_14 = zeros(N_G-1,N_G+1);   % pressure part
    
    % C1 matrix
    C1_gas_11 = zeros(N_G-1,N_G+1);	% radial velocity part
    C1_gas_12 = zeros(N_G-1,N_G+1);   % azimuthal velocity part
    C1_gas_13 = 1i*D0_gas;          % axial velocity part
    C1_gas_14 = zeros(N_G-1,N_G+1);   % pressure part
    
    % C2 matrix
    C2_gas_11 = zeros(N_G-1,N_G+1);	% radial velocity part
    C2_gas_12 = zeros(N_G-1,N_G+1);   % azimuthal velocity part
    C2_gas_13 = zeros(N_G-1,N_G+1);   % axial velocity part
    C2_gas_14 = zeros(N_G-1,N_G+1);   % pressure part
    
    %% radial momentum
    
    % C0 matrix
    C0_gas_21 = (Vr/Dr)*(4/(Rel*(lr-1)^2))*D2_gas + (Vr/Dr)*(4./(Rel*(lr-1)*(rG*(lr-1)+lr+1))).*D1_gas +...
        (-((Vr/Dr)*((4*m_theta^2)./(Rel*(rG*(lr-1)+lr+1).^2))) -(Vr/Dr)*(4./(Rel*(rG*(lr-1)+lr+1).^2)) + 1i*omega).*D0_gas;	% radial velocity part
    C0_gas_22 = -(Vr/Dr)*(1/Rel)*(1i*8*m_theta./((rG*(lr-1)+lr+1).^2)).*D0_gas; % azimuthal velocity part
    C0_gas_23 = zeros(N_G-1,N_G+1);       % axial velocity part
    C0_gas_24 = -2/(Dr*(lr-1))*D1_gas;  % pressure part
    
    % C1 matrix
    C1_gas_21 = -1i*Ug.*D0_gas;     % radial velocity part
    C1_gas_22 = zeros(N_G-1,N_G+1);   % azimuthal velocity part
    C1_gas_23 = zeros(N_G-1,N_G+1);   % axial velocity part
    C1_gas_24 = zeros(N_G-1,N_G+1);   % pressure part
    
    % C2 matrix
    C2_gas_21 = -(Vr/Dr)*(1/Rel)*D0_gas;	% radial velocity part
    C2_gas_22 = zeros(N_G-1,N_G+1);   % azimuthal velocity part
    C2_gas_23 = zeros(N_G-1,N_G+1);   % axial velocity part
    C2_gas_24 = zeros(N_G-1,N_G+1);   % pressure part
    
    %% azimuthal momentum
    
    % C0 matrix
    C0_gas_31 = (Vr/Dr)*(1/Rel)*(8*1i*m_theta./((rG*(lr-1)+lr+1).^2)).*D0_gas;  % radial velocity part
    C0_gas_32= (Vr/Dr)*(4/(Rel*(lr-1)^2))*D2_gas + (Vr/Dr)*(4./(Rel*(lr-1)*(rG*(lr-1)+lr+1))).*D1_gas +...
        (-((Vr/Dr)*((4*m_theta^2)./(Rel*(rG*(lr-1)+lr+1).^2))) - (Vr/Dr)*(4./(Rel*(rG*(lr-1)+lr+1).^2)) + 1i*omega).*D0_gas;    % azimuthal velocity part
    C0_gas_33 = zeros(N_G-1,N_G+1);   % axial velocity part
    C0_gas_34 = -(2*1i*m_theta./(Dr*(rG*(lr-1)+lr+1))).*D0_gas; % pressure part
    
    % C1 matrix
    C1_gas_31 = zeros(N_G-1,N_G+1);   % radial velocity part
    C1_gas_32= -1i*Ug.*D0_gas;      % azimuthal velocity part
    C1_gas_33 = zeros(N_G-1,N_G+1);   % axial velocity part
    C1_gas_34 = zeros(N_G-1,N_G+1);   % pressure part
    
    % C2 matrix
    C2_gas_31 = zeros(N_G-1,N_G+1);   % radial velocity part
    C2_gas_32= -(Vr/Dr)*(1/Rel)*D0_gas; % azimuthal velocity part
    C2_gas_33 = zeros(N_G-1,N_G+1);   % axial velocity part
    C2_gas_34 = zeros(N_G-1,N_G+1);   % pressure part
    
    %% axial momentum
    
    % C0 matrix
    C0_gas_41 = -(2/(lr-1))*Ugp.*D0_gas;    % radial velocity part
    C0_gas_42 = zeros(N_G-1,N_G+1);           % azimuthal velocity part
    C0_gas_43 = (Vr/Dr)*(4/(Rel*(lr-1)^2))*D2_gas + (Vr/Dr)*(4./(Rel*(lr-1)*(rG*(lr-1)+lr+1))).*D1_gas +...
        (-((Vr/Dr)*((4*m_theta^2)./(Rel*(rG*(lr-1)+lr+1).^2))) + 1i*omega).*D0_gas; % axial velocity part
    C0_gas_44 = zeros(N_G-1,N_G+1);   % pressure part
    
    % C1 matrix
    C1_gas_41 = zeros(N_G-1,N_G+1);   % radial velocity part
    C1_gas_42 = zeros(N_G-1,N_G+1);   % azimuthal velocity part
    C1_gas_43 = -1i*Ug.*D0_gas;     % axial velocity part
    C1_gas_44 = -(1i/Dr)*D0_gas;    % pressure part
    
    % C2 matrix
    C2_gas_41 = zeros(N_G-1,N_G+1);   % radial velocity part
    C2_gas_42 = zeros(N_G-1,N_G+1);   % azimuthal velocity part
    C2_gas_43 = -(Vr/Dr)*(1/Rel)*D0_gas;    % axial velocity part
    C2_gas_44 = zeros(N_G-1,N_G+1);   % pressure part
    
    %% Construction of the C0,C1 and C1 matrices for gas
    
    % zero rows at the beginning is to fill the gas boundary and interface
    % conditions.
    
    C0_gas = [zeros(7,4*N_G+3); C0_gas_21 C0_gas_22 C0_gas_23 C0_gas_24(:,1:end-1);...
        C0_gas_31 C0_gas_32 C0_gas_33 C0_gas_34(:,1:end-1);...
        C0_gas_41 C0_gas_42 C0_gas_43 C0_gas_44(:,1:end-1);...
        C0_gas_11 C0_gas_12 C0_gas_13 C0_gas_14(:,1:end-1)];
    
    C1_gas = [zeros(7,4*N_G+3); C1_gas_21 C1_gas_22 C1_gas_23 C1_gas_24(:,1:end-1);...
        C1_gas_31 C1_gas_32 C1_gas_33 C1_gas_34(:,1:end-1);...
        C1_gas_41 C1_gas_42 C1_gas_43 C1_gas_44(:,1:end-1);...
        C1_gas_11 C1_gas_12 C1_gas_13 C1_gas_14(:,1:end-1)];
  
    C2_gas = [zeros(7,4*N_G+3); C2_gas_21 C2_gas_22 C2_gas_23 C2_gas_24(:,1:end-1);...
        C2_gas_31 C2_gas_32 C2_gas_33 C2_gas_34(:,1:end-1);...
        C2_gas_41 C2_gas_42 C2_gas_43 C2_gas_44(:,1:end-1);...
        C2_gas_11 C2_gas_12 C2_gas_13 C2_gas_14(:,1:end-1)];
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Boundary condition at liquid centerline %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % m_theta=0: Axisymmetric mode 
    % m_theta=1: Asymmetric mode     
    
    % radial velocity

    if (m_theta==1)
        tempC0 = [D1_liqCenter zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
        tempC1 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
        tempC2 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    elseif (m_theta==0)
        tempC0 = [D0_liqCenter zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
        tempC1 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
        tempC2 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    end
    
    C0_liq(4*N_L-3,:) = tempC0;
    C1_liq(4*N_L-3,:) = tempC1;
    C2_liq(4*N_L-3,:) = tempC2;
    
    
    % azimuthal velocity
    
    if (m_theta==1)
        tempC0 = [-1i*D0_liqCenter D0_liqCenter zeros(1,N_L+1) zeros(1,N_L)];
        tempC1 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
        tempC2 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    elseif (m_theta==0)
        tempC0 = [zeros(1,N_L+1) D0_liqCenter zeros(1,N_L+1) zeros(1,N_L)];
        tempC1 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
        tempC2 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    end
    
    C0_liq(4*N_L-2,:) = tempC0;
    C1_liq(4*N_L-2,:) = tempC1;
    C2_liq(4*N_L-2,:) = tempC2;
    
    
    % axial velocity
    
    if (m_theta==1)
        tempC0 = [zeros(1,N_L+1) zeros(1,N_L+1) D0_liqCenter zeros(1,N_L)];
        tempC1 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
        tempC2 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    elseif (m_theta==0)
        tempC0 = [zeros(1,N_L+1) zeros(1,N_L+1) D1_liqCenter zeros(1,N_L)];
        tempC1 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
        tempC2 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    end
   
    C0_liq(4*N_L-1,:) = tempC0;
    C1_liq(4*N_L-1,:) = tempC1;
    C2_liq(4*N_L-1,:) = tempC2;
    
    
    % pressure
    
    if (m_theta==1)
        tempC0 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) D0_liqCenter(:,1:end-1)];
        tempC1 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
        tempC2 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    elseif (m_theta==0)
        tempC0 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) D1_liqCenter(:,1:end-1)];
        tempC1 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
        tempC2 = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    end
   
    C0_liq(4*N_L,:) = tempC0;
    C1_liq(4*N_L,:) = tempC1;
    C2_liq(4*N_L,:) = tempC2;
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Boundary condition at gas boundary %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % radial velocity

    tempC0 = [D0_gasBoundary zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];
    tempC1 = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];
    tempC2 = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];
    
    C0_gas(7,:) = tempC0;
    C1_gas(7,:) = tempC1;
    C2_gas(7,:) = tempC2;
    
    % azimuthal velocity
    
    tempC0 = [zeros(1,N_G+1) D0_gasBoundary zeros(1,N_G+1) zeros(1,N_G)];
    tempC1 = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];
    tempC2 = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)]; 
    
   
    C0_gas(6,:) = tempC0;
    C1_gas(6,:) = tempC1;
    C2_gas(6,:) = tempC2;
    
    % axial velocity
    
    tempC0 = [zeros(1,N_G+1) zeros(1,N_G+1) D1_gasBoundary zeros(1,N_G)];
    tempC1 = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];
    tempC2 = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];
    
   
    C0_gas(5,:) = tempC0;
    C1_gas(5,:) = tempC1;
    C2_gas(5,:) = tempC2;
    
    % pressure
    
    tempC0 = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) D1_gasBoundary(:,1:end-1)];
    tempC1 = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];
    tempC2 = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];
   
    C0_gas(4,:) = tempC0;
    C1_gas(4,:) = tempC1;
    C2_gas(4,:) = tempC2;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Combined matrices for the two phases %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % This is excluding the interface conditions. Interface conditions will
    % be added in the next section
    
    C0 = [C0_liq zeros(4*N_L+3,4*N_G+3); zeros(4*N_G+3,4*N_L+3) C0_gas];
    C1 = [C1_liq zeros(4*N_L+3,4*N_G+3); zeros(4*N_G+3,4*N_L+3) C1_gas];
    C2 = [C2_liq zeros(4*N_L+3,4*N_G+3); zeros(4*N_G+3,4*N_L+3) C2_gas];
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Interface conditions %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    %% Jump in Normal stress condition
    
    tempC0_liq = [-4/Rel*D1_liqInterface zeros(1,N_L+1) (1/2/Wel)*(1/(Ulp_int-(1/(lr-1))*Ugp_int))*(-1+m_theta^2)*D0_liqInterface D0_liqInterface(:,1:end-1)];
    tempC1_liq = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    tempC2_liq = [zeros(1,N_L+1) zeros(1,N_L+1) (1/2/Wel)*(1/(Ulp_int-(1/(lr-1))*Ugp_int))*D0_liqInterface zeros(1,N_L)];
    
    tempC0_gas = [4*Vr/(Rel*(lr-1))*D1_gasInterface zeros(1,N_G+1) -(1/2/Wel)*(1/(Ulp_int-(1/(lr-1))*Ugp_int))*(-1+m_theta^2)*D0_gasInterface -D0_gasInterface(:,1:end-1)];
    tempC1_gas = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];
    tempC2_gas = [zeros(1,N_G+1) zeros(1,N_G+1) -(1/2/Wel)*(1/(Ulp_int-(1/(lr-1))*Ugp_int))*D0_gasInterface zeros(1,N_G)];

    C0(4*N_L+1,:) = [tempC0_liq tempC0_gas];
    C1(4*N_L+1,:) = [tempC1_liq tempC1_gas];
    C2(4*N_L+1,:) = [tempC2_liq tempC2_gas];
    
    
    %% Continuity of Shear stress condition-1
    
    tempC0_liq = [-4*Ulpp_int*D0_liqInterface zeros(1,N_L+1) 2i*omega*D1_liqInterface zeros(1,N_L)];
    tempC1_liq = [-omega*D0_liqInterface zeros(1,N_L+1) -2i*Uint*D1_liqInterface zeros(1,N_L)];
    tempC2_liq = [Uint*D0_liqInterface zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    
    tempC0_gas = [(4*Vr/((lr-1)^2))*Ugpp_int*D0_gasInterface zeros(1,N_G+1) -2i*Vr*(omega/(lr-1))*D1_gasInterface zeros(1,N_G)];
    tempC1_gas = [omega*Vr*D0_gasInterface zeros(1,N_G+1) 2i*Vr*(Uint/(lr-1))*D1_gasInterface zeros(1,N_G)];
    tempC2_gas = [-Vr*Uint*D0_gasInterface zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];

    C0(4*N_L+2,:) = [tempC0_liq tempC0_gas];
    C1(4*N_L+2,:) = [tempC1_liq tempC1_gas];
    C2(4*N_L+2,:) = [tempC2_liq tempC2_gas];

    
    %% Continuity of Shear stress condition-2
    
    tempC0_liq = [-1i*m_theta/(rL_interface+1)*D0_liqInterface -D1_liqInterface+(1/(rL_interface+1))*D0_liqInterface zeros(1,N_L+1) zeros(1,N_L)];
    tempC1_liq = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    tempC2_liq = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    
    tempC0_gas = [1i*m_theta*(Vr/(rG_gasInterface*(lr-1)+lr+1))*D0_gasInterface (Vr/(lr-1))*D1_gasInterface-(Vr/(rG_gasInterface*(lr-1)+lr+1))*D0_gasInterface zeros(1,N_G+1) zeros(1,N_G)];
    tempC1_gas = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];
    tempC2_gas = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];

    C0(4*N_L+3,:) = [tempC0_liq tempC0_gas];
    C1(4*N_L+3,:) = [tempC1_liq tempC1_gas];
    C2(4*N_L+3,:) = [tempC2_liq tempC2_gas];
    
    
    %% Continuity of radial velocity condition
    
    tempC0_liq = [-D0_liqInterface zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    tempC1_liq = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    tempC2_liq = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    
    tempC0_gas = [D0_gasInterface zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];
    tempC1_gas = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];
    tempC2_gas = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];

    C0(4*N_L+4,:) = [tempC0_liq tempC0_gas];
    C1(4*N_L+4,:) = [tempC1_liq tempC1_gas];
    C2(4*N_L+4,:) = [tempC2_liq tempC2_gas];
    
    
    %% Continuity of azimuthal velocity condition
    
    tempC0_liq = [zeros(1,N_L+1) -D0_liqInterface zeros(1,N_L+1) zeros(1,N_L)];
    tempC1_liq = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    tempC2_liq = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    
    tempC0_gas = [zeros(1,N_G+1) D0_gasInterface zeros(1,N_G+1) zeros(1,N_G)];
    tempC1_gas = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];
    tempC2_gas = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];

    C0(4*N_L+5,:) = [tempC0_liq tempC0_gas];
    C1(4*N_L+5,:) = [tempC1_liq tempC1_gas];
    C2(4*N_L+5,:) = [tempC2_liq tempC2_gas];
    
    
    %% Continuity of axial velocity condition
    
    tempC0_liq = [-2*Ulp_int*D0_liqInterface zeros(1,N_L+1) 1i*omega*D0_liqInterface zeros(1,N_L)];
    tempC1_liq = [zeros(1,N_L+1) zeros(1,N_L+1) -1i*Uint*D0_liqInterface zeros(1,N_L)];
    tempC2_liq = [zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L+1) zeros(1,N_L)];
    
    tempC0_gas = [(2/(lr-1))*Ugp_int*D0_gasInterface zeros(1,N_G+1) -1i*omega*D0_gasInterface zeros(1,N_G)];
    tempC1_gas = [zeros(1,N_G+1) zeros(1,N_G+1) 1i*Uint*D0_gasInterface zeros(1,N_G)];
    tempC2_gas = [zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G+1) zeros(1,N_G)];

    C0(4*N_L+6,:) = [tempC0_liq tempC0_gas];
    C1(4*N_L+6,:) = [tempC1_liq tempC1_gas];
    C2(4*N_L+6,:) = [tempC2_liq tempC2_gas];
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Solving for the eigenvalues %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    % The non-linear eigenvalue system k^2*C2.a+k*C1.a+C0.a=0 is linearized
    % to the form A.f=k*B.f. Here, f=[k*a; a] where k is the eigenvalue and
    % f is the eigenvector. 'a' vector is the chebyshev coefficients for
    % the radial component of the perturbation quantities
    
    A=[-C1 -C0; eye(size(C1,1),size(C1,2)) zeros(size(C1,1),size(C1,2))];
    B=[C2 zeros(size(C2,1),size(C2,2)); zeros(size(C2,1),size(C2,2)) eye(size(C2,1),size(C2,2))];
    
    [F,k] = eig(B,A);   % eigenvalues being solved
    k=1./diag(k);
    
    ind = find((abs(k-omega))<(omega+1e-6));
    setOfk1 = k(ind);
    setOfF1 = F(:,ind);
    
    if(i==1)
        setOfk = setOfk1;
        setOfF = setOfF1;
    else
        ind1 = find(real(setOfk1)>kFinal(i-1,1));
        setOfk = setOfk1(ind1);
        setOfF = setOfF1(:,ind1);
    end
    
    [M,I] = min(imag(setOfk));  % most unstable (most negative of the complex part of k) eigenvalue is selected (M is the value and I is the index)
    kFinal(i,1) = setOfk(I,1);
    Cheb_coeff(:,i) = setOfF(:,I);
    Wavenumber(i,1) = real(kFinal(i,1));    % real part of the selected eigenvalue is the wavenumber
    GrowthRate(i,1) = -imag(kFinal(i,1));   % negative part of the complex part of the eigenvalur is the growth rate
   
    Error(i,1)=norm(A*Cheb_coeff(:,i)-kFinal(i,1)*B*Cheb_coeff(:,i)); % Residual Error (norm(A.f-k*B.f))
    disp(['Input frequency = ' num2str(omega)]);
    disp(['Growth rate = ' num2str(GrowthRate(i,1))]);
    disp(['Wavenumber = ' num2str(Wavenumber(i,1))]);
    disp(['Residual Error = ' num2str(Error(i,1))]);
end

end