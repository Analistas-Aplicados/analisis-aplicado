%========================================================================
% Curso Análisis Aplicado
% Marzo 2015
% Carlos Dioney Blanco González 
%
% Este programa lleva a cabo el método de gradiente conjugado
% computacional para obtener la solución a un s.e.l. de la forma:
%
%         A*x = b 
%
% El cual es equivalente al problema  de minimización de funciones
% cuadráticas convexas de la siguiente forma
%
%    min  Phi(x) = .5 x'Ax - b'x       A s.p.d.
%   
% Para el caso de minimizar una función tomamos A y b como la matriz
% hessiana de la función y el gradiente respectivamente.
%========================================================================

function x = gradconj2(H, g, tol, maxiter)

%========================================================================
% INPUT:
%  H        - La matriz Hessiana de la función f
%  g        - El gradiente de la función f
% tol       - Tolerancia para la norma del residuo r = A*x - b
%             y para evitar curvatura negativa
% maxiter   - Máximo de iteraciones del método GC
%
% OUTPUT:
%    x      - Solución al sistema de ecuaciones lineales 
%========================================================================


n = length(g);      % Tamaño del vector gradiente
x = zeros(n,1);     % Aproximación inicial
r = H*x - g;        % Obtenemos el residuo
d = -r;             % Nuestra primera dirección en GC
k = 0;              % Contador

% Implementamos el método del gradiente mientras
% nuestro residuo siga siendo mayor a la tolerancia
% o hasta que se alcance el máximo de iteraciones

normr = norm(r);        % Norma del residuo

while ( normr > tol ) && (k < maxiter)
    
    aux    =  H * d;    
    daux   =  d'*aux; 

    if daux <= tol     % Evitamos curvatura negativa
        return
    end

    alpha  =  ( r'*r )/( daux );
    x  	   =  x + alpha*d;
    rvieja =  r;
    r 	   =  r + alpha*aux;
    beta   =  (r'*r) / ( rvieja'*rvieja );
    dvieja =  d;
    d 	   =  -r + beta*dvieja;
    k 	   =  k + 1;
    normr  = norm(r);
    
end