%========================================================================
% Este programa lleva a cabo el método
% BFGS con regiones de confianza
% para obtener una dirección de descenso.
% Curso Análisis Aplicado
% Marzo 2015
% Carlos Dioney Blanco González 
%========================================================================

function [ x, i ] = BFGS_TR( fun, x )

%========================================================================
% INPUT:
%  fun      - La matriz Hessiana de la función f
%  x        - El gradiente de la función f
%
% OUTPUT:
%  p       - Dirección de descenso 
%  i       - Número total de iteraciones
%========================================================================


tol      = 1e-8;             % tolerancia
n        = length(x);        % tamaño del vector
maxiter  = 500;              % máximo de iteraciones
i        = 0;                % iniciamos nuestro contador
[~,g,~]  = feval(fun, x);    % gradiente en el punto inicial
H        = eye(n);           % creamos la primera aproximación de H

deltamax = 1000;             % Para regiones de confianza

while (norm(g) > tol && i < maxiter) 
    
       p     = CG_TR( H, g, deltamax, tol); % Obtenemos una dirección con BFGS
             
       xn    = x + p;               % Nuevo punto 
    [~,g2,~] = feval(fun, xn);      % Nuevo gradiente
       s     = xn - x;              
       y     = g2 - g;
       
    % Construimos la nueva aproximación de H con BFGS
    H = H + ( (y*y') / (y'*s) ) - ( (H*(s*s')*H) / (s'*H*s) );
    x = xn;          % Actualizamos   
    g = g2;
    i = i + 1;
end

