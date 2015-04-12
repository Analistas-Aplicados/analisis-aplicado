%========================================================================
% Este programa lleva a cabo el m�todo
% BFGS con regiones de confianza
% para obtener una direcci�n de descenso.
% Curso An�lisis Aplicado
% Marzo 2015
% Carlos Dioney Blanco Gonz�lez 
%========================================================================

function [ x, i ] = BFGS_TR( fun, x )

%========================================================================
% INPUT:
%  fun      - La matriz Hessiana de la funci�n f
%  x        - El gradiente de la funci�n f
%
% OUTPUT:
%  p       - Direcci�n de descenso 
%  i       - N�mero total de iteraciones
%========================================================================


tol      = 1e-8;             % tolerancia
n        = length(x);        % tama�o del vector
maxiter  = 500;              % m�ximo de iteraciones
i        = 0;                % iniciamos nuestro contador
[~,g,~]  = feval(fun, x);    % gradiente en el punto inicial
H        = eye(n);           % creamos la primera aproximaci�n de H

deltamax = 1000;             % Para regiones de confianza

while (norm(g) > tol && i < maxiter) 
    
       p     = CG_TR( H, g, deltamax, tol); % Obtenemos una direcci�n con BFGS
             
       xn    = x + p;               % Nuevo punto 
    [~,g2,~] = feval(fun, xn);      % Nuevo gradiente
       s     = xn - x;              
       y     = g2 - g;
       
    % Construimos la nueva aproximaci�n de H con BFGS
    H = H + ( (y*y') / (y'*s) ) - ( (H*(s*s')*H) / (s'*H*s) );
    x = xn;          % Actualizamos   
    g = g2;
    i = i + 1;
end

