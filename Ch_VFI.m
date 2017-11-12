%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binghamton University                         %
% ECON634 Advanced Macroeconomics               %
% Luis Chancí (lchanci1@binghamton.edu)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program does the Value Function Iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[v0,g,idx]=Ch_VFI(u,v0,beta,PI,a)
  num_a = length(a);
  e1 = 1;
    while e1 >.0001; %1e-06
          v = u + beta * repmat(permute((PI*v0),[3 2 1]),[num_a 1 1]);
          [vfn,idx] = max(v,[],2);        
          e1  = max(max(abs(permute(vfn,[3 1 2]) - v0)));
          v0  = permute(vfn,[3 1 2]);
    end
    idx = permute(idx,[3 1 2]);
    g   = a(idx);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%