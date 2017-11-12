%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binghamton University                         %
% ECON634 Advanced Macroeconomics               %
% Luis Chancí (lchanci1@binghamton.edu)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script does the VFI and the linear       %
% interpolation method                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u            = Ch_U(r,w,a0,z,sigma);
v0           = zeros(m, num_a0);
V            = zeros(m,num_a);

[v0,g0,idx0] = Ch_VFI(u,v0,beta,PI,a0);

for i=1:m
	V(i,:)   = interp1(a0, v0(i,:), a);
end
u            = Ch_U(r,w,a,z,sigma);
v            = u + beta * repmat(permute(PI*V,[3 2 1]),[num_a 1 1]);
[vfn,idx]    = max(v,[],2); 
v0           = permute(vfn,[3 1 2]);
idx          = permute(idx,[3 1 2]);
g            = a(idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
