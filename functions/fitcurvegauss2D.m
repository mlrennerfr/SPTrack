function [pestimates, fval,exitflag,sse]  = fitcurvegauss2D(p, data)

% Call fminsearch with a random starting point.
[maxX, maxY]=size(data);
%start_point=[p(6),p(7)];
%posx=p(1);
%posy=p(2);
%intens=p(4);
%offset=p(5);

options=optimset('Display','off');

[pestimates, fval,exitflag] = fminsearch(@gaussfun, p,options);

sse = gaussfun(pestimates);

% gaussfun accepts curve parameters as inputs and outputs sse,
% the sum of squares error 

    function sse = gaussfun(p)
        
        %FittedCurve = A * exp( -x / (2*sigma^2) );
        for y=1:maxY
            for x=1:maxX
               FittedCurve(x,y)=p(4)*exp(-2*((x-p(1))^2/(p(6)^2))-2*((y-p(2))^2/(p(7)^2)))+p(5);
               Error(x,y) = FittedCurve(x,y) - data(x,y);
            end
        end
        
        sse = sum(sum(Error .^ 2));

    end

end


%f=h*exp(-2*((x-xo)^2/(widthx^2))-2*((y-yo)^2/(widthy^2)))+b
