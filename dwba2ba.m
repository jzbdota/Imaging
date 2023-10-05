function F = dwba2ba(G,tranA)

% regroupG
[fs1,fs2,fs3] = size(G);
Gr = nan(fs1,fs2,fs3*2,class(G));
Grr = nan(fs1,fs2,fs3*4,class(G));
Gtmp = nan(fs1+1,fs2+1,fs3,class(G));
if isa(G,'gpuArray')
    Gr = single(Gr);
    Grr = single(Grr);
    Gtmp = single(Gtmp);
end
Gtmp(1:fs1,1:fs2,:) = G;
Gtmp(end,1:fs2,:) = G(1,:,:);
Gtmp(1:fs1,end,:) = G(:,1,:);
Gtmp(end,end,:) = G(1,1,:);
Gtmp = rot90(Gtmp,2);
Gtmp(end,:,:) = [];
Gtmp(:,end,:) = [];
Gr(:,:,1:fs3) = G;
Gr(:,:,(fs3+1):end) = Gtmp;
Grr(:,:,1:(fs3*2)) = real(Gr);
Grr(:,:,(fs3*2+1):end) = imag(Gr);
Grr = permute(Grr,[3,1,2]);
Grr = reshape(Grr,fs3*4,fs1*fs2);

Fr = tranA*Grr;

% degroupF
fs3d = size(Fr,1);
fs3 = fs3d/2;
Fr = reshape(Fr,fs3d,fs1,fs2);
Fr = permute(Fr,[2,3,1]);
F = Fr(:,:,1:fs3) + 1i*Fr(:,:,(fs3+1):end);
return