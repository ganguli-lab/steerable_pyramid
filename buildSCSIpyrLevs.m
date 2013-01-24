% [PYR, INDICES] = buildSCSIpyrLevs(LODFT, LOGRAD, XRCOS, YRCOS, ANGLE, HEIGHT, NBANDS)
%
% Recursive function for constructing levels of a steerable pyramid.  This
% is called by buildSCSIpyr, and is not usually called directly.

function [pyr,pind] = buildSCSIpyrLevs(lodft,log_rad,Xrcos,Yrcos,angleMap,ht,nbands)

% if the height is 0, we are at the end of the pyramid
if (ht <= 0)

		% return the low-passed DFT in the image domain
		lo0 = ifft2(ifftshift(lodft));
		pyr = real(lo0(:));
		pind = size(lo0);

else

		% bands contains the coefficients for the different orientations at this frequency
		bands = zeros(numel(lodft), nbands);

		% indices for these bands
		bind = zeros(nbands,2);

		% shift origin of LUT by 1 octave.
		Xrcos = Xrcos - log2(2);

		% LUT
		lutsize = 1024;
		Xcosn = pi*(-(2*lutsize+1):(lutsize+1))/lutsize;
		order = nbands-1;

		%% divide by sqrt(sum_(n=0)^(N-1)  cos(pi*n/N)^(2(N-1)) )
		%% Thanks to Patrick Teo for writing this out :)
		const = (2^(2*order))*(factorial(order)^2)/(nbands*factorial(2*order));

		% Ycosn = sqrt(const) * (cos(Xcosn)).^order;
		% analytic version: only take one lobe
		alfa=	mod(pi+Xcosn,2*pi)-pi;
		Ycosn = 2*sqrt(const) * (cos(Xcosn).^order) .* (abs(alfa)<pi/2);

		% build the high freq mask...
		himask = pointOp(log_rad, Yrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);

		for b = 1:nbands

				% angle mask for oriented filter
				anglemask = pointOp(angleMap, Ycosn, Xcosn(1)+pi*(b-1)/nbands, Xcosn(2)-Xcosn(1));

				% compute the 2D FFT masked by the high-freq mask and the angled (oriented) mask
				banddft = ((-1i)^(nbands-1)) .* lodft .* anglemask .* himask;

				% store these coefficients
				band = ifft2(ifftshift(banddft));
				bands(:,b)=band(:);

				% store the size of the band
				bind(b,:)  = size(band);

		end

		% generate a new low-freq mask (one octave reduced)
		YIrcos = abs(sqrt(1.0 - Yrcos.^2));
		lomask = pointOp(log_rad, YIrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);

		% filter the 2D FFT with the lomask
		lodft = lomask .* lodft;

		% recursively build the rest of the pyramid
		[npyr,nind] = buildSCSIpyrLevs(lodft, log_rad, Xrcos, Yrcos, angleMap, ht-1, nbands);

		% return the pyramid and the indices
		pyr = [bands(:); npyr];
		pind = [bind; nind];

end
