function varargout = edgen(varargin)
% generates electron density matrix and support matrix (3D)
% assuming sample has the same dimension in the x,y directions

if ~nargin
    szxy = 256;
    szz = 64;
else
    szxy = varargin{1};
    szz = varargin{2};
end
gr = importrandobj(szxy,szz);
support = mask3d(szxy,szz);
gr(~support) = 0;

% st = dbstack;
% disp(st.name)

if nargout
    varargout{1}=gr;
    varargout{2}=support;
end

function gr = importrandobj(ftsize,thickness)

x = -ftsize/2:ftsize/2-1;
x = x+.5;
[X,Y,Z] = ndgrid(x,x,x);

defz = 8;

sizerotcube = 8; % size of rot cube is 8*2
% locrotcube = [-16,16,0]; % z symmetric
locrotcube = [-16,16,-8 + defz]; % sitting on the substrate

sizecube = 10;
% loccube = [26,-16,0];
loccube = [26,-16,-5 + defz];

sizerect1 = [4,26,6];
% locrect1 = [21,11,0];
locrect1 = [21,11,-3 + defz];

sizerect2 = [16,4,6];
% locrect2 = [-21,-26,0];
locrect2 = [-21,-26,-3 + defz];

sizerect3 = [8,8,24];
% locrect3 = [0,0,0];
locrect3 = [0,0,-12 + defz];

sizesphere = 8; % radius
% locsphere = [-16,0,0];
locsphere =[-16,0,-8 + defz];

%   sphere
dis = (X-locsphere(1)).^2+(Y-locsphere(2)).^2+(Z-locsphere(3)).^2;
objsphere = dis<sizesphere^2;

%   rect1 rect2 and rect3
rect1 = rect(X-locrect1(1),sizerect1(1)).*rect(Y-locrect1(2),sizerect1(2)).*rect(Z-locrect1(3),sizerect1(3));
rect2 = rect(X-locrect2(1),sizerect2(1)).*rect(Y-locrect2(2),sizerect2(2)).*rect(Z-locrect2(3),sizerect2(3));
rect3 = rect(X-locrect3(1),sizerect3(1)).*rect(Y-locrect3(2),sizerect3(2)).*rect(Z-locrect3(3),sizerect3(3));
objrect1 = logical(rect1);
objrect2 = logical(rect2);
objrect3 = logical(rect3);

%   cube
cube = rect(X-loccube(1),sizecube).*rect(Y-loccube(2),sizecube).*rect(Z-loccube(3),sizecube);
objcube = logical(cube);

%   rotated cube
Xt = X(:)';Yt = Y(:)';Zt = Z(:)';
TMP = [Xt;Yt;Zt];
%%%%%%%
TMPr = rotz(45)*TMP;
% TMPr = rotx(45)*TMPr;
% TMPr = rotx(45)*TMP;
%%%%%%%
Xr = TMPr(1,:);Yr = TMPr(2,:);Zr = TMPr(3,:);
Xr = reshape(Xr',ftsize,ftsize,ftsize);
Yr = reshape(Yr',ftsize,ftsize,ftsize);
Zr = reshape(Zr',ftsize,ftsize,ftsize);
Xr = Xr - locrotcube(1);
Yr = Yr - locrotcube(2);
Zr = Zr - locrotcube(3);

m1 = Xr < sizerotcube & Xr > -sizerotcube;
m2 = Yr < sizerotcube & Yr > -sizerotcube;
m3 = Zr < sizerotcube & Zr > -sizerotcube;
objrotcube = m1.*m2.*m3;

% gr = objsphere.*2 + objrect1 + objrect2 + objrect3 + objcube + objrotcube;
gr = objsphere + objrect1 + objrect2 + objrect3 + objcube + objrotcube;
gr = gr(:,:,(ftsize/2-thickness/2+1):(ftsize/2+thickness/2));
return

function result = mask3d(ftsize,thickness)
% ftsize and thickness should be multiple of 4

if mod(thickness,4) || mod(ftsize,4)
    error('ftsize or thickness is not multiple of 4');
end

thickhalf = thickness/2;
fthalf = ftsize/2;
    
realsize = fthalf;
x = -fthalf:fthalf-1;
x = x+.5;
z = -thickhalf:thickhalf-1;
z = z +.5;
[X,Y,Z] = ndgrid(x,x,z);

%    ellipsoid mask
dis = sqrt(X.^2+Y.^2+(ftsize/thickness)^2*Z.^2);
result = dis<realsize/2;

% rectangular mask
% result = rect(X,fthalf).*rect(Y,fthalf).*rect(Z,thickhalf);
mask2 = rect(X-15,4).*rect(Y+6,4); % defect in the mask, necessary
mask2 = ~mask2;
result = result.*mask2;

result = logical(result);

function y = rect(x,D)
% rectangle function
if nargin == 1
    D = 1;
end
x = abs(x);
y = double(x<D/2);
y(x==D/2) = 0.5;
return

