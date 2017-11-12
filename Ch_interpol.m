function [V]= Ch_interpol(x,m,v0,a,opt)
V = [];
if opt==1     % Linear

elseif opt==2 % Cubic Splines
    for i=1:m
        Polynm = spline(a, v0(i,:));
        V=[V;ppval(Polynm,x)];
    end
% elseif opt==2
end
