%========================================================================
% Este programa lleva a cabo el m�todo
% BFGS con b�squeda lineal
% para obtener la soluci�n al sistema de ecuaciones 
% A*x = b
% Curso An�lisis Aplicado
% Marzo 2015
% Carlos Dioney Blanco Gonz�lez 
%========================================================================

function [ p, i ] = BFGS_lineal( fun, x )

%========================================================================
% INPUT:
%  fun      - La funci�n a minimizar
%  x        - El punto inicial
%
% OUTPUT:
%  p       - Soluci�n del problema
%  i       - N�mero total de iteraciones
%========================================================================


tol      = 1e-8;             % tolerancia
maxiter  = 50;               % m�ximo de iteraciones
i        = 0;                % iniciamos nuestro contador
[fx,g,~] = feval(fun, x);    % funci�n y gradiente en el punto inicial
n        = length(x);              
B        = eye(n);           % creamos la primera aproximaci�n de H


while (norm(g) > tol && i < maxiter)
    
    p    = -B \ g;           % Obtenemos nuestra direcci�n
    
    % El tama�o �ptimo del paso con b�squeda lineal
    alfa = linesch_sw( x, fx, g, p, fun, 0.0001, 0.9, 2); 
    
       xn    = x + alfa*p;          % Nuevo punto 
    [~,g2,~] = feval(fun, xn);      % Nuevo gradiente
       s     = xn - x;              
       y     = g2 - g;
       
    % Construimos la nueva aproximaci�n de H con BFGS
    B = B + ( (y*y') / (y'*s) ) - ( (B*(s*s')*B) / (s'*B*s) );
    x = xn;          % Actualizamos   
    g = g2;
    i = i + 1;
end

p = x;  % Regresamos la soluci�n

