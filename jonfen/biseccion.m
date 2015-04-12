function [ alfa ] = biseccion (x, f, g, p, funder,c1,c2)
%--------------------------------------------------------------------------

falla = 0;                 %  1 si el procedimiento hizo mas de 20 pasos
numfg = 0;                 %  numero de evaluaciones de f y g
alfa  = 1;  gTp = g'*p;    %  derivada direccional
%
% c1 = 0.0001; c2 = 0.9;
while numfg < 20
    xp = x + alfa*p;
    [ fp, gp, ~ ] = feval( funder, xp );
    numfg = numfg + 1;
    if fp <= f + c1*alfa*gTp && abs(gp'*p) <= c2*abs(gTp)
        x = xp; f = fp; g = gp;
        return
    else
        alfa = alfa/2;
    end
end
falla = 1;
