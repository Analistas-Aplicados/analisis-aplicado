%========================================================================
% Curso An�lisis Aplicado
% Marzo 2015
% Carlos Dioney Blanco Gonz�lez 
%
% Este programa lleva a cabo el m�todo de gradiente conjugado
% computacional para obtener la soluci�n a un s.e.l. de la forma:
%
%         A*x = b 
%
% El cual es equivalente al problema  de minimizaci�n de funciones
% cuadr�ticas convexas de la siguiente forma
%
%    min  Phi(x) = .5 x'Ax - b'x       A s.p.d.
%   
% Para el caso de minimizar una funci�n tomamos A y b como la matriz
% hessiana de la funci�n y el gradiente respectivamente.
%========================================================================

function x = gradconj2(H, g, tol, maxiter)

%========================================================================
% INPUT:
%  H        - La matriz Hessiana de la funci�n f
%  g        - El gradiente de la funci�n f
% tol       - Tolerancia para la norma del residuo r = A*x - b
%             y para evitar curvatura negativa
% maxiter   - M�ximo de iteraciones del m�todo GC
%
% OUTPUT:
%    x      - Soluci�n al sistema de ecuaciones lineales 
%========================================================================


n = length(g);      % Tama�o del vector gradiente
x = zeros(n,1);     % Aproximaci�n inicial
r = H*x - g;        % Obtenemos el residuo
d = -r;             % Nuestra primera direcci�n en GC
k = 0;              % Contador

% Implementamos el m�todo del gradiente mientras
% nuestro residuo siga siendo mayor a la tolerancia
% o hasta que se alcance el m�ximo de iteraciones

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