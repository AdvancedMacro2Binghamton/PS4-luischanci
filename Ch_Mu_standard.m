%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binghamton University                         %
% ECON634 Advanced Macroeconomics               %
% Luis Chancí (lchanci1@binghamton.edu)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script iterates over the distribution    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 e2 = 1; 
 while e2>=1e-06;
        [emp_ind, a_ind, mass] = find(Mu);
        MuNew  = zeros(size(Mu));
        for ii = 1:length(emp_ind)
            apr_ind = idx(emp_ind(ii), a_ind(ii));   
            MuNew (:, apr_ind) = MuNew (:, apr_ind) + ... 
                  (PI(emp_ind(ii), :) * mass(ii))';
        end
        e2 = max(max(abs(MuNew  - Mu)));
        Mu = MuNew;
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%