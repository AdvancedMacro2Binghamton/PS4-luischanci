%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binghamton University                         %
% PhD in Economics                              %
% ECON634 Advanced Macroeconomics               %
% Fall 2017                                     %
% Luis Chancí (lchanci1@binghamton.edu)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% OPTIONS (METHOD):
% 0 = VALUE FUNCTION ITERATION
% 1 = POLICY FUNCTION ITERATION
% 2 = GRID SEARCH
% 3 = INTERPOLATION - LINEAR
% 4 = INTERPOLATION - CUBIC SPLINE
   
  Option = 0 ;
  
% PARAMETERS
beta    = 0.99;     sigma   = 2.0;    
alpha   = 1/3;      delta   = 0.025;
rho     = 0.5;      sigma_e = 0.2; 
tol     = 0.05;

% DISCRETIZATION:
m                  = 5;
[z, PI]            = TAUCHEN(m, rho, sigma_e, 3);
z                  = exp(z'); 
[eig_vec, eig_val] = eig(PI');
PI_ss              = bsxfun(@rdivide,eig_vec,sum(eig_vec));

% ASSETS, CAPITAL AND LABOR:
a_lo   = 0;            
a_hi   = 80;           
num_a  = 500;   % If want to see resuts for Opt=4 on the same day, use only 50. 
num_a0 = 0.1*num_a;   % Initial Grid for interpolation
a      = linspace(a_lo, a_hi, num_a);
a0     = linspace(a_lo, a_hi, num_a0);
K_min  = 20; K_max = 40;
L      = (z*PI_ss(:,1)); 

e0 = 1; itera = 1; t0 = tic;
while (e0 >tol)        % (e0 >0.01)||(itera<=1000)
 K_guess = (K_min + K_max) / 2;
 r       = alpha*(K_guess ^ (alpha - 1))*(L ^ (1 - alpha))+(1-delta);
 w       = (1-alpha)*(K_guess ^ alpha)*(L^(-alpha));
 
 % I. VALUE FUNCTION (SOLUTION):
 
 if     Option == 0
        u          = Ch_U(r,w,a,z,sigma);
        v0         = zeros(m, num_a);
        [v0,g,idx] = Ch_VFI(u,v0,beta,PI,a);
 elseif Option == 1
        Ch_PFI
 elseif Option == 2
        Ch_Grid_Search
 elseif Option == 3
        Ch_Interpol_Linear
 elseif Option == 4
        % (i) Initial guess (short grid).
        u            = Ch_U(r,w,a0,z,sigma);
        v0           = zeros(m, num_a0);
        [V0,G0,idx0] = Ch_VFI(u,v0,beta,PI,a0);
        c0           = bsxfun(@plus, r*a0',permute(z,[1 3 2])*w);
        [V1,G1,VP]   = Ch_Interpol_Spline(V0,G0,idx0,c0,a0,PI,sigma,beta);
                   
        % (ii) Interpolation for the complete grid
        u            = Ch_U(r,w,a,z,sigma);
        v_guess      = ppval(VP,a);
        [V2,G2,idx2] = Ch_VFI(u,v_guess,beta,PI,a);
        c            = bsxfun(@plus, r*a',permute(z,[1 3 2])*w);
        [v0,g,vp]    = Ch_Interpol_Spline(V2,G2,idx2,c,a,PI,sigma,beta);      
 end
 
% II. ITERATION OVER THE DISTRIBUTION
 Mu = ones(size(g))/numel(g);
 if Option <= 3
     Ch_Mu_standard
 else
     Ch_Mu_Spline
 end

% III. MARKETS CLEAR
 aggsav = sum(sum(Mu.*g));
 if (aggsav - K_guess)>0
    K_min = K_guess; elseif (aggsav - K_guess)<0; K_max = K_guess; 
 end
 e0 = abs(aggsav - K_guess); itera = itera + 1;
 display (['K = ', num2str(K_guess)])
 display (['Aggregate desired wealth = ', num2str(aggsav)]);
end;  t = toc(t0);           % Now we can plot :)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%