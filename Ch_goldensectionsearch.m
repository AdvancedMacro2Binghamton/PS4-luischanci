function [val, x_opt] = goldensectionsearch( f, a, b )
% Copyright (c) 2009, Katarzyna Zarnowiec, all rights reserved
f = @(x) -f(x);
epsilon=0.000001;               % accuracy value
iter= 50;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k=0;                            % number of iterations
x1=a+(1-tau)*(b-a);             % computing x values
x2=a+tau*(b-a);
f_x1=f(x1);                     % computing values in x points
f_x2=f(x2);
while ((abs(b-a)>epsilon) && (k<iter))
    k=k+1;
    if(f_x1<f_x2)
        b=x2;
        x2=x1;
        x1=a+(1-tau)*(b-a);
        f_x1=f(x1);
        f_x2=f(x2);
    else
        a=x1;
        x1=x2;
        x2=a+tau*(b-a);
        f_x1=f(x1);
        f_x2=f(x2);
    end
    k=k+1;
end
if(f_x1<f_x2)
    x_opt = x1;
    val   = -f_x1;
else
    x_opt = x2;
    val   = -f_x2;
end
end