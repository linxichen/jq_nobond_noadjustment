function retval = chebyeval(a,x)
% usage: y = chebyeval(a,x)
% description:  given input Chebyshev vector a and vector
% x of absicssas, this routine returns the vector of values of 
% corresponding Chebyshev polynomial at the coordinates of x.
 
% local variables: 
% j,k:  index variables
% m: size of coefficient array
% n: size of node array
% z: 2*x(j)
% xpar,xtmp: temporaries

m = length(a);
n = length(x);
retval = zeros(size(x));
for j = 1:n % loop over entries of x 
  % use Clenshaw's method to evaluate the polynomial
  z = 2*x(j);
  xpar = zeros(1,3);
  for k = 1:m
    xtmp = z*xpar(3) - xpar(2) + a(k);
    xpar = [xpar(2:3), xtmp];
  end;
  retval(j) = (xpar(3)-xpar(1))/2;
end;
