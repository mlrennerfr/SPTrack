function [a,da,chi] = fitgaussLM(options,a,y,dy)
%function [a,da,chi] = fitgaussLM(options,a,y,dy)
% Non-linear least-square fit of gauss function using the Levenberg/Marquard algorithm.
%	a()     - fitting parameter
%	y()     - data array
%	chi   - reduced chi-squared
%	a()   - fitted parameter
%   da()  - error in the fitting parameters
% modified from marqogaussian.m
% ts <6.1994> WJ & GAB <8.2000>
% MR mar 09 for SPTrack v4.0
% --------------------------------------------------------------------

% check input parameter and prepare internal variables
npar = prod(size(a));
[ymax,xmax] = size(y);
npnt = xmax*ymax;
delta     = 0.01;
lambda    = 0.0001;
chimin    = options(1)*(npnt-npar);
dchimin   = options(2)*(npnt-npar);
dpmin     = options(3);
maxtry    = options(4)*npar;
maxlambda = options(5);
dpar      = 1;
oldchi    = 1e18;

if nargin==3
   dy = ones(ymax,xmax);
   issigma=0;
else
   dy = 1 ./dy;
   issigma=1;
end

dyvec = dy(:);
par=a;
deriv=zeros(npnt,npar);
if issigma
	b=(fgauss(par,xmax,ymax)-y).*dy;
else
	b=(fgauss(par,xmax,ymax)-y);
end
b=b(:); 
chi = b'*b;
tries=1; 
docontinue=1;
   
while docontinue

    % calculate function and chi squared
	ytry = fgauss(par,xmax,ymax);   
	if issigma
      b=(ytry-y).*dy; 
	else
      b=(ytry-y);  
    end
    b=b(:);
    newchi = b'*b;
   
	%actualize parameter if chi sqared decreased 
	if newchi<chi 
   	    oldchi = chi;
		chi = newchi;
		lambda = 0.1*lambda;
		a = par;
	else
	   lambda = 10*lambda;
	end

	%check for break condition
	if (tries>=maxtry) | (chi<chimin) | (lambda>maxlambda) | (min(a./dpar')>1/dpmin) | (oldchi-chi<dchimin)
      lambda = 0;
      docontinue = 0;
  	else
      tries=tries+1;
    end
     
	%calculate the derivatives
	if chi==newchi
   	   for i = 1:npar
	      par = a;
	      if a(i) == 0
   			 par(i) = delta;
             h = (fgauss(par, xmax, ymax)-ytry) / delta;
          else
			 par(i) = a(i)*(1+delta);
			 h = (fgauss(par, xmax, ymax)-ytry) / (delta*a(i));
          end
          h = h(:);
          if issigma
             deriv(:,i) = h .* dyvec;
          else
             deriv(:,i) = h;
          end
	   end
	end

	%calculate alpha and beta matrices
    if issigma
      h = (y-ytry) .* dy;
    else
      h = (y-ytry);
    end
	h = h(:);
	alpha = deriv' * deriv;
	alpha = alpha + lambda*diag(diag(alpha));
	beta  = deriv' * h;
   
    if lambda ~= 0
		% solve linear equation and actualize parameter
		if rank(alpha)==length(alpha)
			dpar = alpha \ beta;
			par = a + dpar';
        else
   		    docontinue=0;
      end
   end
end

chi = newchi/(npnt-npar);
if rank(alpha)==length(alpha)
  alpha = inv(alpha);
  da    = (sqrt(chi) * sqrt(diag(alpha)))';
else
  da = a;
end

% eof

