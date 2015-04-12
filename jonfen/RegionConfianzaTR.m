%========================================================================
% Regiones de Confianza 
% Curso An�lisis Aplicado
% Marzo 2015
% Carlos Dioney Blanco Gonz�lez
%========================================================================

function [ x, f, i] = ...
RegionConfianzaTR( x, H, deltamax, delta, etha, tol, maxiter)

%========================================================================
% INPUT: 
%    x 	   - El punto inicial
% deltamax - El tama�o m�ximo de la regi�n de confianza
% delta0   - Regi�n de confianza inicial,  delta0 est� en ( 0, deltamax)
%   etha   - Un valor entre ( 0, 0.25 ]
%   tol	   - Tolerancia del m�todo
%  maxiter - M�ximo de iteraciones
%
% OUTPUT:
%    x	   - Valor donde se minimiza la funci�n
%    i     - N�mero de iteraciones necesarias
%    f     - Valor de la funci�n evaluada en el punto x
%========================================================================

i = 0;  	% Iniciamos nuestro contador

% Definimos nuestro modelo, buscamos una direcci�n TR

m = @(p, g, B, f) (.5 * p' * B * p) + (g'* p) + f;

%========================================================================
%   m  -  Representa nuestro modelo, donde
%
%   p  -  Direcci�n de descenso obtenida con GC
%   g  -  El gradiente de la funci�n evaluado en x
%   B  -  La Hessiana de la funci�n evaluado en x
%   f  -  Funci�n evaluada en x
%========================================================================

[~, g, ~]   = rosenbrock(x);		% Obtenemos nuestro gradiente
normg       = norm(g);              % Evaluamos su norma 
cero        = zeros(length(g),1);	% Creamos un vector de ceros de n*1

% Iniciamos el ciclo mientras la norma del gradiente sea mayor a la
% tolerancia y no lleguemos al m�ximo de iteraciones.

while (normg > tol) && (i < maxiter)
    
    % Obtenemos p_k con gradiente conjugado con delta0, la Hessiana de f y
    % gradiente de f  ( f es nuestra funci�n rosenbrock )

    [f, g, ~] = rosenbrock(x);    
    
    % evaluamos el radio ( rho )
    rho1 = rosenbrock(x) - rosenbrock( x + p );
    rho2 = m( cero, g, H, f) - m( p, g, H, f);
    rho  = rho1 / rho2;
   
    % checamos las condiciones para rho y ver c�mo movemos delta 
    if  rho < .25       
        delta = .25 * norm(p);
        
    elseif (rho > .75) && (norm(p) == delta)
        delta = min( 2*delta, deltamax);
    end
    
    if rho > etha
        x = x - p;
    end
    i = i + 1;
    normg = norm(g);
end