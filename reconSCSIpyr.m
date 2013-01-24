% RES = reconSCSIpyr(PYR, INDICES, LEVS, BANDS, TWIDTH)
%
% The inverse of buildSCSIpyr: Reconstruct image from its complex steerable pyramid representation,
% in the Fourier domain.

function res = reconSCSIpyr(pyr, indices, Nsc, Nor, levs, bands, twidth)

%%------------------------------------------------------------
%% DEFAULTS:

if ~exist('levs'),
  levs = 'all';
end

if ~exist('bands')
  bands = 'all';
end

if ~exist('twidth'),
  twidth = 1;
elseif (twidth <= 0)
  fprintf(1,'Warning: TWIDTH must be positive.  Setting to 1.\n');
  twidth = 1;
end

%%------------------------------------------------------------

pind = indices;
%Nsc = log2(pind(1,1)/pind(end,1)); % incorrect!!
%Nor = (size(pind,1)-2)/Nsc;

for nsc = 1:Nsc,
    firstBnum = (nsc-1)*Nor+2;

%% Re-create analytic subbands
    dims = pind(firstBnum,:);
    ctr = ceil((dims+0.5)/2);
    ang = mkAngle(dims, 0, ctr);
    ang(ctr(1),ctr(2)) = -pi/2;
    for nor = 1:Nor,
      nband = (nsc-1)*Nor+nor+1;
      ind = pyrBandIndices(pind,nband);
      ch = pyrBand(pyr, pind, nband);
      ang0 = pi*(nor-1)/Nor;
      xang = mod(ang-ang0+pi, 2*pi) - pi;
      amask = 2*(abs(xang) < pi/2) + (abs(xang) == pi/2);
      amask(ctr(1),ctr(2)) = 1;
      amask(:,1) = 1;
      amask(1,:) = 1;
      amask = fftshift(amask);
      ch = ifft2(amask.*fft2(ch));    % "Analytic" version
      %f = 1.000008;  % With this factor the reconstruction SNR goes up around 6 dB!
      f = 1;
      ch = f*0.5*real(ch); % real part
      pyr(ind) = ch;
    end     % nor
end         % nsc

res = reconSFpyr(pyr, indices, levs, bands, twidth);
