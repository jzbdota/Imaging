function varargout=disperr(A,B,DISPLAY)
% display the error between A and B
% A and B can be matrices with any same dimension (1D, 2D, 3D,...)
% A and B can be real and complex numbers
% see also: norm

narginchk(2,3);
nargoutchk(0,1);
if nargin<3
    DISPLAY = true;
end
% A = round(A,4); B = round(B,4);
err = norm(A(:)-B(:))/norm(A(:)+B(:));
if DISPLAY
    disp(err);
end
if nargout
    varargout{1}=err;
end
end