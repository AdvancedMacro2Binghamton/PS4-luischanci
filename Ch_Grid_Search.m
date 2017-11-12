%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binghamton University                         %
% ECON634 Advanced Macroeconomics               %
% Luis Chancí (lchanci1@binghamton.edu)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script conducts a grid search            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 u       = Ch_U(r,w,a,z,sigma);
 v0      = zeros(m, num_a);
 e1 = 1;
   while e1 > tol
   	v         = u + beta * repmat(permute(PI*v0,[3 2 1]),[num_a 1 1]);
    [vfn,idx] = max(v,[],2); 
    e1        = max(max(abs(permute(vfn,[3 1 2]) - v0)));
    v0        = permute(vfn,[3 1 2]);
   end;
 idx     = permute(idx,[3 1 2]);
 g       = a(idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%