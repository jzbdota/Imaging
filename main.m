function varargout = main
% 10/14/2023
% Yi Yang, yiy042@ucsd.edu
% HIO-ER
% Including Gaussian noise and Beamstop
% See also: edgen.m, slscattering.m

nargoutchk(0,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controling flags
HIOONLY = false;
BEAMSTOP = true;
NOISE = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
if NOISE
    SNR = 10;
end

disp(datetime);
shrink_threshold = .13; % shrink wrap method threshold, [0,1]
totaliter = 1000; % total number of iterations
mfname = mfilename; % name of the m file

global ftsize thickness;
ftsize = 256; thickness = 64;

if ~HIOONLY
    [gro,support] = edgen(ftsize,thickness);
    % Add gaussian profile, details are given in the function mygauss
    gaussian_profile = mygauss();
    gro = gro.*gaussian_profile;
    
    if exist('SNR','var')
        [Im_sqrt,A,tranA] = slscattering(gro,SNR);
    else
        [Im_sqrt,A,tranA] = slscattering(gro);
    end
    save(mfname,'gro','A','tranA','support','Im_sqrt')
    if BEAMSTOP
        beamstopind = blockpixels(1);
        save(mfname,'beamstopind','-append');
    end

else
    load(mfname,'gro','A','tranA','support','Im_sqrt');
    if BEAMSTOP
        load(mfname,'beamstopind');
    end
end

INDICATOR = false;
R = nan(1,totaliter,'single'); % Error initialization

gr = ifftn(dwba2ba(Im_sqrt.*exp(1i*2*pi*rand(size(Im_sqrt))),tranA));

% GPU accelaration
R = gpuArray(R);
gr = gpuArray(single(gr));
supportgpu = gpuArray(support);
Im_sqrt = gpuArray(single(Im_sqrt));
A = gpuArray(single(A));
tranA = gpuArray(single(tranA));
if BEAMSTOP
    beamstopind = gpuArray(beamstopind);
end

tic
for itern = 1:totaliter
    Btmp = fftn(gr);
    grn = gr;
    Gtmp = ba2dwba(Btmp,A);
    if BEAMSTOP
        Im_sqrt(beamstopind) = abs(Gtmp(beamstopind));
    end
    gr = ifftn(dwba2ba(Im_sqrt.*exp(1i*angle(Gtmp)),tranA));
    gr(~supportgpu) = grn(~supportgpu) - .9*gr(~supportgpu);
    R(itern) = disperr(fftn(gr),fftn(grn),false);

    if INDICATOR
        INDICATOR = false;
        gr = real(gr);
        gr(gr<0)=0;
        gr(~supportgpu) = 0;
    end
    if ~mod(itern,60)
        gr = abs(gr);
        gr(~supportgpu) = 0;
        INDICATOR = true;
    end

    % Shrinkwrap hybrid
    if ~mod(itern,100)
        sigma= -itern*2/totaliter + 3.3;
        sigma = gather(sigma);
        grsmooth = imgaussfilt3(abs(gr),sigma);% Bug if sigma is gpuArray
        supportgpu = grsmooth > shrink_threshold;

    % for better shrinkwrap
        supportgpu(:,:,1:(ftsize/32-1)) = 0; % depends on mask3d
        supportgpu(:,:,(thickness-ftsize/32+2):end) = 0;
    end

    if ~mod(itern,totaliter/10) && size(dbstack,1) == 1
        name = ['Iteration ' num2str(itern)];
        disp(name);
    end
end
toc

gr(~supportgpu) = 0;
gr = abs(gr);
gr = gather(gr);
save(mfname,'gr','-append');
R = gather(R);
save(mfname,'R','-append');

Rf = rferr(gr,gro);
disp(['Rf is : ',num2str(Rf)]);

if size(dbstack,1) == 1
    figure;semilogy(R);
    xlabel('Iteration#');
    ylabel('Convergence');
    legend(strcat('Rf is ',num2str(round(Rf*100,2)),'%'));
    ISO = .4;
    gr = flip(gr,3);
    gro = flip(gro,3);
    figure('visible','on');isosurface  (gro,ISO); title('Original');
    figure; isosurface  (gr,ISO); title('Reconstructed');
end

outcell = {gr,Rf,R};
if nargout
    varargout = cell(nargout,1);
    for ii = 1:nargout
        varargout{ii}=outcell{ii};
    end
end
return

function gp = mygauss()
% Gaussian Profile in the y direction
% As footprint in the x direction is much bigger than sample size
% We assume the beam in the x and z direction is uniform
% fuzzy logic toolbox is required.

global ftsize thickness;

x = 1:ftsize;
x = x - .5;

gp_vector = gaussmf(x,[32,ftsize/2]);

y_support = 96:160;
gp_support = gp_vector(y_support);
mean_gp = mean(gp_support);

gp_vector = gp_vector / mean_gp;
gp = repmat(gp_vector,[ftsize,1,thickness]);

% figure;plot(y_support,gp_support/mean_gp);

return

function err = rferr(A,B)
% abs(fftn(A)), abs(fftn(B))
err = disperr(abs(fftn(A)),abs(fftn(B)),false);
return

function beamstopind = blockpixels(Dhalf)
    
global ftsize thickness;
beamstopind = [];
pixvecxy = [1:Dhalf+1, (ftsize-Dhalf+1):ftsize];
pixvecz = [1:Dhalf+1, (thickness-Dhalf+1):thickness];
for ii = pixvecxy
    for jj = pixvecxy
        for kk = pixvecz
            ind = sub2ind([ftsize ftsize thickness],ii,jj,kk);
            beamstopind = [beamstopind,ind];
        end
    end
end
return
