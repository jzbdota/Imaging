function G = ba2dwba(F,A)
% regroupF
[fs1,fs2,fs3] = size(F);
Fr = nan(fs1,fs2,fs3*2,class(F));
if isa(F,'gpuArray')
    Fr = single(Fr);
end
Fr(:,:,1:fs3) = real(F);
Fr(:,:,(fs3+1):end) = imag(F);

Fr = permute(Fr,[3,1,2]);
Fr = reshape(Fr,fs3*2,fs1*fs2);

Gr = A*Fr;
% degroupG
[gs3,~] = size(Gr);
lssize = gs3/4;
Gr = reshape(Gr,gs3,fs1,fs2);
Gr = permute(Gr,[2,3,1]);
Grr = Gr(:,:,1:(lssize*2)) + 1i*Gr(:,:,(lssize*2+1):end);
G = Grr(:,:,1:lssize);
return