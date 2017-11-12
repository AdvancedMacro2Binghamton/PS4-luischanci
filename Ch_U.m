%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binghamton University                         %
% ECON634 Advanced Macroeconomics               %
% Luis Chancí (lchanci1@binghamton.edu)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is for the utility function      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = Ch_U(r,w,a,z,s)
 c      = bsxfun(@plus, bsxfun(@minus,r*a',a),permute(z,[1 3 2])*w);
 u      = (c.^(1-s))./(1-s); 
 u(c<0) = -Inf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%