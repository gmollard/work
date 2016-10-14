function [xvect, resvect, nit] = newton( fun, dfun, x0, tol, nmax )
% NEWTON Find a zero of a nonlinear scalar function.
%   [XVECT] = NEWTON(FUN,DFUN,X0,TOL,NMAX) finds a zero of the differentiable 
%   function FUN using the Newton method and returns a vector XVECT containing 
%   the successive approximations of the zero (iterates). DFUN is the derivative of FUN.
%   FUN and DFUN accept real scalar input x and return a real scalar value; 
%   FUN and DFUN can also be inline objects. X0 is the initial guess.
%   TOL is the tolerance on error allowed and NMAX the maximum number of iterations.
%   The stopping criterion based on the difference of successive iterates is used.
%   If the search fails a warning message is displayed.
%   
%   [XVECT,RESVECT,NIT] = NEWTON(FUN,DFUN,X0,TOL,NMAX) also returns the vector
%   RESVECT of residual evaluations for each iterate, and NIT the number of iterations.
%   Note: the length of the vectors is equal to ( NIT + 1 ).
%

xvect = [ x0 ];
res = tol+1;
resvect = [res];
i=0;
nit = [i];
x= x0;
while (i<=nmax) && (abs(fun(xvect(i+1))) > tol)
    i=i+1;
    nit = [nit i];
    x = x - (fun(x)/dfun(x));
    xvect = [ xvect x ];
    resvect = [ resvect xvect(i+1)-xvect(i) ];
    res = resvect(i+1);
end
    

return