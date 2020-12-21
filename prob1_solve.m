% [a, b] = solvef(2);
a = 1;
b = 0;
x = 0:5;
x = x.*pi./10; 
a*x-b

% a = 2;
% b = 0*pi/10*a;
% f(a,b)

A = 1/6*[sum(x.^2), sum(x); -sum(x), 6]
[u,s,v] = svd(A);
s
ss = 2/(s(1)+s(2))

function [a,b] = solvef(i)
    x = i:5;
    x = x.*pi./10;
    gx = 1 - cos(x);
    x2 = x.^2;
    temp = sum(gx)/(6-i);
    a_n = sum(x.*gx) - sum(x)*temp;
    a_d = sum(x2) - 1/(6-i) * sum(x)^2;
    a = a_n / a_d
    b = (a * sum(x) - sum(gx)) / (6-i)
    
    res = f(a,b)
    
    if a * x(1) - b < 0
        fprintf("no sol\n")
    end
    
    if i > 0 && a * (i-1)*pi/10 - b > 0
        fprintf("no sol\n")
    end
    
end


function [ga, gb] = gf(a, b)
    x = 0:5;
    x = x.*pi./10;
    ga = 0;
    gb = 0;
    for i = 1:6
        rel = relu(a*x(i)-b);
        reg = rel - g(x(i));
        grel = grelu(a*x(i)-b);
        ga = ga + x(i)*grel*reg;
        gb = gb + grel*reg;
    end
    ga = ga / 6;
    gb = - gb / 6;
end

function res = f(a, b)
    x = 0:5;
    x = x.*pi./10;
    res = 0;
    for i = 1:6
       re = relu(a*x(i)-b) - g(x(i));
       res = res + re^2;
    end
    res = res / 12;
end

function re = relu(x)
    re = max(x,0);
end

function gre = grelu(x)
    if x > 0
        gre = 1;
    else
        gre = 0;
    end
end

function gx = g(x)
    gx = 1 - cos(x);
end