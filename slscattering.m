function [result,A,tranA] = slscattering(gr,varargin)
% single layer DWBA scattering
% assuming gr has the same dimension in the x,y directions
% assuming gr is NOT mirror symmetric regarding to z direction
% Incident angle should be bigger than outgoing
% if incident is smaller, swap incident and outgoing

% See also: edgen

if nargin == 2
    SNR = varargin{1};
end

[ftsize,~,thickness]=size(gr);
fthalf = ftsize/2;
thickhalf = thickness/2;
[qxmesh,qymesh] = ndgrid(0:ftsize-1);

ki = 0:thickhalf-1;
angles = .01 + .02*ki;
% angles = angles * 4; % bigger angles for higher z resolution

lssize = 0;  % number of angle pairs
anglein = nan(thickhalf^2,1); % incident angle array
angleout = nan(thickhalf^2,1); % outgoing array
for ii = ki
    for ff = ki
        if (ii+ff)>=thickhalf || ii<ff % ii<ff necessary for FriedelD(A3)
            continue
        end
        lssize = lssize + 1;
        anglein(lssize) = angles(ii+1);
        angleout(lssize) = angles(ff+1);
    end
end
anglein((lssize+1):end) = [];
angleout((lssize+1):end) = [];
k0 = .5/sind(angles(1)); % incident wavevector magnitude

kzi_vac = k0*sind(anglein); % incident kz in vacuum
kzi_sub = k0*sqrt(cukaxray.sin^2-cosd(anglein).^2); % incident kz in substrate
kzf_vac = k0*sind(angleout); % outgoing kz in vacuum
kzf_sub = k0*sqrt(cukaxray.sin^2-cosd(angleout).^2);
Ri = (kzi_vac-kzi_sub)./(kzi_vac+kzi_sub);
Rf = (kzf_vac-kzf_sub)./(kzf_vac+kzf_sub);
D1 = ones(size(Ri),class(Ri));
D2 = Rf;
D3 = Ri;
D4 = Ri.*Rf;
qz1 = kzf_vac + kzi_vac;
qz2 = -kzf_vac + kzi_vac;
qz3 = -qz2;
qz4 = -qz1;

rowA = 1:lssize;
Asize1 = thickness;
A1 = sparse(rowA,mod(round(qz1),Asize1)+1,D1,lssize,Asize1);
A2 = sparse(rowA,mod(round(qz2),Asize1)+1,D2,lssize,Asize1);
A3 = sparse(rowA,mod(round(qz3),Asize1)+1,D3,lssize,Asize1);
A4 = sparse(rowA,mod(round(qz4),Asize1)+1,D4,lssize,Asize1);
A1 = single(full(A1)); A2 = single(full(A2));
A3 = single(full(A3)); A4 = single(full(A4));
A = newmatrix(A1,A2,A3,A4);
tranA = (A'*A)\A';

%%%%%%%%%%% The angetheta resolution and rotation resolution is GOOD
angletheta = linspace(-.1,.1,ftsize*2);

anglesup1 = linspace(-.02,0,ftsize/8);
anglesup2 = -anglesup1;
angletheta = [angletheta,anglesup1,anglesup2];

angletheta = unique(angletheta);
qy0 = sind(angletheta);
qpnorm = (sqrt(2)+.1)*fthalf/max( abs(qy0) );

result = nan(ftsize,ftsize,lssize);
tic
for lsind = 1:lssize
    if ~mod(lsind,40)
        disp(lsind);
        toc
        tic
    end

    qx0 = cosd(angleout(lsind)).*cosd(angletheta) - cosd(anglein(lsind));
    qy0 = sind(angletheta).*cosd(angleout(lsind));
    qx0 = qpnorm * qx0; qy0 = qpnorm * qy0;
    qxccd = []; qyccd = [];

    for jj = 1:length(qx0)
        anglerot = linspace(0,180,180); % not good, sometimes good
        qxtmp = qx0(jj)*cosd(anglerot)-qy0(jj)*sind(anglerot);
        qytmp = qx0(jj)*sind(anglerot)+qy0(jj)*cosd(anglerot);
        qxccd = [qxccd,qxtmp];
        qyccd = [qyccd,qytmp];
    end
    % delete redundant qx qy
    qxccdp = qxccd';
    qyccdp = qyccd';
    qccdmesh = [qxccdp,qyccdp];
    qccdmesh = unique(qccdmesh,'rows');
    qxccd = squeeze(qccdmesh(:,1));
    qyccd = squeeze(qccdmesh(:,2));
    qxccd = qxccd';
    qyccd = qyccd';
    % delete big qx qy
    tempind = find(qxccd < -fthalf-2 | qxccd > fthalf+1);
    qxccd(tempind) = []; qyccd(tempind) = [];
    tempind = find(qyccd < -fthalf-2 | qyccd > fthalf+1);
    qxccd(tempind) = []; qyccd(tempind) = [];

    [ftccd1,ftccd2,ftccd3,ftccd4] = myft3dwbagpu(gr,qxccd,qyccd, ...
    qz1(lsind),qz2(lsind),qz3(lsind),qz4(lsind));
    ftccd_dwba = D1(lsind)*ftccd1 + D2(lsind)*ftccd2 ...
    + D3(lsind)*ftccd3 + D4(lsind)*ftccd4;
% no diffuse scattering from rough surface
    ftccd_dwba = abs(ftccd_dwba);
% add noise according to SNR. SNR of diffraction pattern, which is square
% of Fourier transform amplitudes
    if exist('SNR','var')
%         test_dwba = ftccd_dwba;
        ftccd_dwba = ftccd_dwba.^2;
        ftccd_dwba = ftccd_dwba .* (1+1/SNR*randn(size(ftccd_dwba),'like',ftccd_dwba));
        ftccd_dwba = sqrt(ftccd_dwba);
    end
%     disperr(test_dwba,ftccd_dwba,true);
%     pause(1000);
    
    ft_dwba = griddata(qxccd,qyccd,ftccd_dwba,qxmesh-fthalf,qymesh-fthalf,'natural');
    ft_dwba = fftshift(ft_dwba);
    result(:,:,lsind) = ft_dwba;
end
toc
return

function Af = FriedelD(A)

Af = A;
Af(:,2:end) = A(:,end:-1:2);
return

function rpMAT = newmatrix(D1,D2,D3,D4)
% detail is given in DWBA_Matrix Summary_ver2.docx

Dplus = D1 + D2 + D3 + D4;
Dminus = FriedelD(Dplus);
rdp = real(Dplus); idp = imag(Dplus);
rdm = real(Dminus); idm = imag(Dminus);
rpMAT = [rdp,-idp;rdm,idm;idp,rdp;idm,-rdm];
% save('nmat.mat',"Dplus","Dminus","rdp","idp","idm","rdm");
% disp('Data saved!');
% pause(1000)
return

function [r1,r2,r3,r4] = myft3dwbagpu(X,qx,qy,qz1,qz2,qz3,qz4)
% qzs are scalars
% maxmeshsize = 32 is the best for RTX 2080/2070(8GB VRAM)

maxmeshsize = gpuArray(single(32*1024*1024));

X = single(X);
qx = single(qx);
qy = single(qy);
qz1 = single(qz1);
qz2 = single(qz2);
qz3 = single(qz3);
qz4 = single(qz4);

N1 = single(size(X,1)); N2 = single(size(X,2)); N3 = single(size(X,3));

[J1,J2,J3] = ndgrid(0:N1-1,0:N2-1,0:N3-1);
% delete 0s to speed up
X=X(:);J1=J1(:);J2=J2(:);J3=J3(:);
tmp = (X==0);
X(tmp)=[];J1(tmp)=[];J2(tmp)=[];J3(tmp)=[];
% disp(size(qx));pause(1000)

Yt1 = gpuArray(nan(size(qx),'single'));
Yt2 = gpuArray(nan(size(qx),'single'));
Yt3 = gpuArray(nan(size(qx),'single'));
Yt4 = gpuArray(nan(size(qx),'single'));
J1gpu = gpuArray(J1(:));
J2gpu = gpuArray(J2(:));
J3gpu = gpuArray(J3(:));
Xgpu = gpuArray(X(:));
N1 = gpuArray(N1); N2 = gpuArray(N2); N3 = gpuArray(N3);

qxgpu = gpuArray(qx);
qygpu = gpuArray(qy);
qz1gpu = gpuArray(qz1);
qz2gpu = gpuArray(qz2);
qz3gpu = gpuArray(qz3);
qz4gpu = gpuArray(qz4);

matqsize = round(maxmeshsize/length(Xgpu));
qlength = length(qxgpu);
qnsteps = ceil(qlength/matqsize);
j3mesh = repmat(J3gpu,1,matqsize);
Xmesh = repmat(Xgpu,1,matqsize);
tempw31 = exp(-1i*2*pi*qz1gpu*j3mesh/N3);
tempw32 = exp(-1i*2*pi*qz2gpu*j3mesh/N3);
tempw33 = exp(-1i*2*pi*qz3gpu*j3mesh/N3);
tempw34 = exp(-1i*2*pi*qz4gpu*j3mesh/N3);
for qstep = 0:qnsteps-1
    if qstep<(qnsteps-1)
        vectemp = (qstep*matqsize+1):(qstep+1)*matqsize;
    else
        vectemp = (qstep*matqsize+1):qlength;
        j3mesh = repmat(J3gpu,1,length(vectemp));
        Xmesh = repmat(Xgpu,1,length(vectemp));
        tempw31 = exp(-1i*2*pi*qz1gpu*j3mesh/N3);
        tempw32 = exp(-1i*2*pi*qz2gpu*j3mesh/N3);
        tempw33 = exp(-1i*2*pi*qz3gpu*j3mesh/N3);
        tempw34 = exp(-1i*2*pi*qz4gpu*j3mesh/N3);
    end
    qxtemp = qxgpu(vectemp);
    qytemp = qygpu(vectemp);
    
    [qxmesh,j1mesh] = meshgrid(qxtemp,J1gpu);
    [qymesh,j2mesh] = meshgrid(qytemp,J2gpu);
    
    tmp = Xmesh.*exp(-1i*2*pi*(qxmesh.*j1mesh/N1+qymesh.*j2mesh/N2));
    fqs = tmp.*tempw31;
    fqs = sum(fqs,1);
    Yt1(vectemp) = squeeze(fqs);
    fqs = tmp.*tempw32;
    fqs = sum(fqs,1);
    Yt2(vectemp) = squeeze(fqs);
    fqs = tmp.*tempw33;
    fqs = sum(fqs,1);
    Yt3(vectemp) = squeeze(fqs);
    fqs = tmp.*tempw34;
    fqs = sum(fqs,1);
    Yt4(vectemp) = squeeze(fqs);
end
r1 = double(gather(Yt1));
r2 = double(gather(Yt2));
r3 = double(gather(Yt3));
r4 = double(gather(Yt4));
return

