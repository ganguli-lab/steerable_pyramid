% builds a steerable, complex, scale invariant pyramid
% Niru Maheswaranathan
% Wed Jan 23 18:58:31 2013

% This is a modified version of buildSCFpyr, that constructs a scale-invariant
% complex-valued steerable pyramid  using Hilbert-transform pairs
% of filters.  Note that the imaginary parts will *not* be steerable.

function [pyr, pind] = buildSCFpyr(im, ht, order, twidth)

%-----------------------------------------------------------------
%% DEFAULTS:

% maximum possible height of the pyramid
max_ht = floor(log2(min(size(im)))) - 2;

% get pyramid height
if exist('ht', 'var') ~= 1
  ht = max_ht;
else
  if (ht > max_ht)
    error('Cannot build pyramid higher than %d levels.',max_ht);
  end
end

% get the order of the pyramid (number of orientations - 1)
if exist('order', 'var') ~= 1
  order = 3;
elseif ((order > 15)  | (order < 0))
  fprintf(1,'Warning: ORDER must be an integer in the range [0,15]. Truncating.\n');
  order = min(max(order,0),15);
else
  order = round(order);
end
nbands = order+1;

% twidth parameter (?)
if exist('twidth', 'var') ~= 1
  twidth = 1;
elseif (twidth <= 0)
  fprintf(1,'Warning: TWIDTH must be positive.  Setting to 1.\n');
  twidth = 1;
end

%-----------------------------------------------------------------

% image dimensions and center
dims = size(im);
ctr = ceil((dims+0.5)/2);

% ramp is a linearly increasing gradient from -1 to 1
[xramp,yramp] = meshgrid( ([1:dims(2)]-ctr(2))./(dims(2)/2), ...
    ([1:dims(1)]-ctr(1))./(dims(1)/2) );

% angle of the ramps
angleMap = atan2(yramp,xramp);

% log_rad is a radial mask
log_rad = sqrt(xramp.^2 + yramp.^2);
log_rad(ctr(1),ctr(2)) =  log_rad(ctr(1),ctr(2)-1);
log_rad  = log2(log_rad);

%% Radial transition function (a raised cosine in log-frequency), used to generate mask:
[Xrcos,Yrcos] = rcosFn(twidth,(-twidth/2),[0 1]);
Yrcos = sqrt(Yrcos);
YIrcos = sqrt(1.0 - Yrcos.^2);

% build the low-frequency mask, keeps low-freq FFT and zeros out high frequencies
lo0mask = pointOp(log_rad, YIrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);

% compute a 2D Fourier transform of the image
imdft = fftshift(fft2(im));

% mask out the high frequencies using our mask
lo0dft =  imdft .* lo0mask;

% compute steerable pyramid coefficients using the 2D low-passed FFT of our image
[pyr,pind] = buildSCSIpyrLevs(lo0dft, log_rad, Xrcos, Yrcos, angleMap, ht, nbands);

% build a high frequency mask (keeps the high frequency residuals)
hi0mask = pointOp(log_rad, Yrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);

% mask out the low frequencies using the high freq. mask
hi0dft =  imdft .* hi0mask;

% convert the high freq. content back to the image domain
hi0 = ifft2(ifftshift(hi0dft));

% return the pyramid coefficients (high freq residuals followed by the rest of the pyramid)
pyr = [real(hi0(:)) ; pyr];

% return the pyramid indices matrix (high freq indices followed by the rest of the pyramid)
pind = [size(hi0); pind];
