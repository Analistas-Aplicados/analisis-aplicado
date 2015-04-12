%========================================================================
% Método de Newton con Gradiente Conjugado
%
% Este programa utiliza el método de Newton con gradiente conjugado para
% resolver el problema de 
%
% min  Phi(x) = .5 x'Ax - b'x       A s.p.d.
%
% Curso Análisis Aplicado
% Marzo 2015
% Carlos Dioney Blanco González
%========================================================================

function [ x, i, norma] = NewtonGC( fun, x, maxiter)

n     = length(x);
tol   = 1.0e-8;
maxCG = 2*n;
i     = 1;

[f1, g, H]  = feval(fun, x);
     norma  = norm(g);
       

while (norma > tol) && (i < maxiter)    
    
    p           = gradconj2( H, -g, tol, maxCG);
    alpha       = 1;
    [f2, g2, ~] = feval( fun, (x + alpha*p));
   
    alpha = biseccion(x, f1, f2, g, g2, p, fun);    
    x     = x + alpha*p;
    
    [f1, g, H] = feval(fun, x);
    i          = i + 1;
    norma      = norm(g);
    
end

end



function alpha =  biseccion(x, f1, f2, g, g2, p, fun)

c1    = 10e-4;
c2    = .9; 
alpha = 1;

    while  ( f2 > f1 + c1*alpha*g'*p ) || ( abs(p'*g2) > c2*abs(p'*g) )
    
        alpha       = alpha/2;
        [f2, g2, ~] = feval(fun, (x + alpha*p));
        
    end
end
