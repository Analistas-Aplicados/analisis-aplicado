% Este programa realiza el m�todo de Newton con 
% gradiente conjugado aplicado a la funci�n rosenbrock
%
% Carlos Dioney Blanco Gonz�lez

function [x,i,norma] = NewtonGC1(fun, x, maxiter)
%definimos nuestros par�metros
n = length(x);
c1 = 10e-4;
c2 = .9;
tol = 1.0e-6;
tol2 = 1.0e-12;
maxCG = 2*n;
[f, g, H] = feval(fun, x);
norma = norm(g);
i = 0;
fprintf(1,'    iter      alfa        norm(grad)         f        \n\n');

while (norma > tol) && (i < maxiter)
    %obtenemos nuestra direcci�n por medio de GC
    p = gradconj( H, -g, tol, maxCG, tol2);
  
   
    %aqu� partimos nuestra alpha si no se cumplen las condiciones de wolfe
    
    [alpha, x, ~, ~, ~, ~] = ...
           linesch_sw(x, f, g, p, fun, c1, c2, 0);
    
      %actualizamos nuestro punto X_km nuestro valor de la funci�n, el
    %gradiente y la Hessiana
   
    [f, g, H] = feval(fun, x);
    i = i + 1;
    norma = norm(g);
    
    fprintf(1,'    %3.0f     %5f     %5.5f        %5.5f    \n',i,alpha,norma,f); 
end
    if (i == maxiter)
    fprintf(1, 'El proceso alcanz� el m�ximo n�mero de iteraciones. \n\n' );
    else
    fprintf(1, 'Se realiz� el proceso con �xito. \n\n');
    end
disp('Se alcanz� el m�nimo en:');    