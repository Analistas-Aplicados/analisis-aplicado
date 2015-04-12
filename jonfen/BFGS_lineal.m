%========================================================================
% Este programa lleva a cabo el método
% BFGS con búsqueda lineal
% para obtener la solución al sistema de ecuaciones 
% A*x = b
% Curso Análisis Aplicado
% Marzo 2015
% Carlos Dioney Blanco González 
%========================================================================

function [ p, i ] = BFGS_lineal( fun, x )

%========================================================================
% INPUT:
%  fun      - La función a minimizar
%  x        - El punto inicial
%
% OUTPUT:
%  p       - Solución del problema
%  i       - Número total de iteraciones
%========================================================================


tol      = 1e-8;             % tolerancia
maxiter  = 50;               % máximo de iteraciones
i        = 0;                % iniciamos nuestro contador
[fx,g,~] = feval(fun, x);    % función y gradiente en el punto inicial
n        = length(x);              
B        = eye(n);           % creamos la primera aproximación de H


while (norm(g) > tol && i < maxiter)
    
    p    = -B \ g;           % Obtenemos nuestra dirección
    
    % El tamaño óptimo del paso con búsqueda lineal
    alfa = linesch_sw( x, fx, g, p, fun, 0.0001, 0.9, 2); 
    
       xn    = x + alfa*p;          % Nuevo punto 
    [~,g2,~] = feval(fun, xn);      % Nuevo gradiente
       s     = xn - x;              
       y     = g2 - g;
       
    % Construimos la nueva aproximación de H con BFGS
    B = B + ( (y*y') / (y'*s) ) - ( (B*(s*s')*B) / (s'*B*s) );
    x = xn;          % Actualizamos   
    g = g2;
    i = i + 1;
end

p = x;  % Regresamos la solución

