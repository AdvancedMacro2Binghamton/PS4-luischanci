%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binghamton University                         %
% ECON634 Advanced Macroeconomics               %
% Luis Chancí (lchanci1@binghamton.edu)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script does the Policy Function Iteration%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 u       = Ch_U(r,w,a,z,sigma);
 v0      = zeros(m, num_a);
 e1 = 1;
   while e1 >1e-06
   	v         = u + beta * repmat(permute(PI*v0,[3 2 1]),[num_a 1 1]);
    [vfn,idx] = max(v,[],2); 
    e1        = max(max(abs(permute(vfn,[3 1 2]) - v0)));
    v0        = permute(vfn,[3 1 2]);
    idx       = permute(idx,[3 1 2]);
    g         = a(idx);
    Q         = makeQmatrix(idx, PI);
    c         = bsxfun(@plus, bsxfun(@minus,r*a,g),z'*w);
    ufn       = (c.^(1-sigma))./(1-sigma); ufn = ufn(:);
    W         = v0(:);
    for     j = 1:30
            V = ufn + beta*Q*W; W = V;
    end
   v0         = reshape(W,m,num_a); 
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
