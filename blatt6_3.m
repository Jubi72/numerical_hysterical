X = [-20 -10 0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 10 20];
n = size(X,2);
eps = 1e-8;
kmax = 50;
Y = zeros([2*n n]);
for i = 1:1
for j = 1:1
    x = [2;1.40];
    k = 0;
    while norm(f(x)) > eps && k < kmax
        Deltax = J(x) \ f(x);
        if Deltax == Inf
            break
        end
        x = x + Deltax;
        Y(2*i-1,j) = x(1);
        Y(2*i  ,j) = x(2);
        k = k + 1;
        if norm(f(x)) < eps
            disp(f(x))
        end
    end
end
end
disp(Y);

function [f]=f(x)
f = [x(1)*x(1) + x(2)*x(2) - 6; x(2) - x(1)*x(1)];
end
    
function [J]=J(x)
if size(x,1) == 2 && size (x,2) == 1
    J = [2*x(1) 2*x(2); 2*x(1) -1];
else
    J = 0;
end
end