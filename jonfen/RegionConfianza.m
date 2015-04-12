%========================================================================
% Regiones de Confianza 
% Curso Análisis Aplicado
% Marzo 2015
% Carlos Dioney Blanco González
%========================================================================

function [ x, f, i] = RegionConfianza( x, deltamax, delta, etha, tol, maxiter)

%========================================================================
% INPUT: 
%    x 	   - El punto inicial
% deltamax - El tamaño máximo de la región de confianza
% delta0   - Región de confianza inicial,  delta0 está en ( 0, deltamax)
%   etha   - Un valor entre ( 0, 0.25 ]
%   tol	   - Tolerancia del método
%  maxiter - Máximo de iteraciones
%
% OUTPUT:
%    x	   - Valor donde se minimiza la función
%    i     - Número de iteraciones necesarias
%    f     - Valor de la función evaluada en el punto x
%========================================================================

i = 0;  	% Iniciamos nuestro contador

% Definimos nuestro modelo, buscamos una dirección con GC.
% y la función de rosenbrock

m = @(p, g, B, f) (.5 * p' * B * p) + (g'* p) + f;

%========================================================================
%   m  -  Representa nuestro modelo, donde
%
%   p  -  Dirección de descenso obtenida con GC
%   g  -  El gradiente de la función evaluado en x
%   B  -  La Hessiana de la función evaluado en x
%   f  -  Función evaluada en x
%========================================================================


tol2        = 10e-8;                % Tolerancia para GC
[~, g, ~]   = rosenbrock(x);		% Obtenemos nuestro gradiente
normg       = norm(g);              % Evaluamos su norma 
cero        = zeros(length(g),1);	% Creamos un vector de ceros de n*1

% Iniciamos el ciclo mientras la norma del gradiente sea mayor a la
% tolerancia y no lleguemos al máximo de iteraciones.

while (normg > tol) && (i < maxiter)
    
    % Obtenemos p_k con gradiente conjugado con delta0, la Hessiana de f y
    % gradiente de f  ( f es nuestra función rosenbrock )

    [f, g, H] = rosenbrock(x);
    p = gradconj2( H, g, tol, 2*length(g), tol2 );
    
    % evaluamos el radio ( rho )
    rho1 = rosenbrock(x) - rosenbrock( x + p );
    rho2 = m( cero, g, H, f) - m( p, g, H, f);
    rho  = rho1 / rho2;
   
    % checamos las condiciones para rho y ver cómo movemos delta 
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
    rosenbrock(x)
end