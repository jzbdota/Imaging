classdef cukaxray
% X-ray constants, Si constants and Au constants
    properties(Constant)
% Below is x ray data
% Cu K alpha
        xraye = 8.04; % kev
        xraylambda = 1.54; % Angstrom
        xrayk0 = 4.078; % inverse Angstrom
       
% Below is Au data
       aucrit = rad2deg( sqrt(2*cukaxray.audelta) ); % degree
       audelta = 49.6e-6;
       aubeta = 5.11e-6;
       aurerho = 131.5; % *10^10cm^-2
       aun = 1 - cukaxray.audelta + 1i*cukaxray.aubeta;
% Below is Si data
       sicrit = rad2deg( sqrt(2*cukaxray.sidelta) ); % degree
       sidelta = 7.56e-6;
       sibeta = 1.73e-7;
       sirerho = 20; % *10^10cm^-2
       sin = 1 - cukaxray.sidelta + 1i*cukaxray.sibeta;

% Below is PS data
        pscrit = rad2deg( sqrt(2*cukaxray.psdelta) );
        psdelta = 3.5e-6;
        psbeta = 4.9e-9;
        psrerho = 9.5;
        psn = 1 - cukaxray.psdelta+1i*cukaxray.psbeta;
    end  
end