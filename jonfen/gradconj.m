% Carlos Dioney Blanco González
% Este programa lleva a cabo el método
% del gradiente conjugado

function x = gradconj(H, g, tol, maxiter)

n = length(g);
x = zeros(n,1);
r = g - H*x;
d = r;
k = 1;
rr=r'*r;

while ( sqrt(rr) > tol && k < maxiter)
    
    aux = H * d; 
    daux = d'*aux;
    alpha = ( rr )/( daux );
    x = x + alpha*d;
    rvieja = r;
    r = r - alpha*aux;
    rr = r'*r;
    beta = (rr) / ( rvieja'*rvieja );
    dvieja = d;
    d = r + beta*dvieja;
    k = k + 1;
        
end